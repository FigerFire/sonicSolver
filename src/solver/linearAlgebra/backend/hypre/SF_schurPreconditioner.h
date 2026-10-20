#pragma once
/// @file SF_schurPreconditioner.h
/// @brief HYPRE FGMRES 的串行 Schur 回调适配；MPI/HYPRE 对象只存在于后端。
#include "SF_hypreControls.h"
#include "solver/linearAlgebra/SF_blockSchur.h"
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <mpi.h>
#include <memory>

namespace SF::LinearAlgebra {
/// @brief 可重复写入的串行 IJ 向量桥接。
class SchurVector {
public:
    explicit SchurVector(int size):indices_(size) {
        for(int i=0;i<size;++i)indices_[i]=i;
        hypreCheck(HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,size-1,&ij_),"Schur vector create");
        hypreCheck(HYPRE_IJVectorSetObjectType(ij_,HYPRE_PARCSR),"Schur vector type");
        write(std::vector<double>(size,0));
    }
    ~SchurVector(){if(ij_)HYPRE_IJVectorDestroy(ij_);}
    SchurVector(const SchurVector&)=delete;
    SchurVector& operator=(const SchurVector&)=delete;
    void write(const std::vector<double>& values) {
        if(values.size()!=indices_.size())throw std::runtime_error("Schur vector size mismatch.");
        hypreCheck(HYPRE_IJVectorInitialize(ij_),"Schur vector initialize");
        hypreCheck(HYPRE_IJVectorSetValues(ij_,(int)indices_.size(),indices_.data(),values.data()),"Schur vector write");
        hypreCheck(HYPRE_IJVectorAssemble(ij_),"Schur vector assemble");
        hypreCheck(HYPRE_IJVectorGetObject(ij_,reinterpret_cast<void**>(&par_)),"Schur vector object");
    }
    std::vector<double> read(HYPRE_ParVector vector) const {
        std::vector<double> result(indices_.size());
        hypreCheck(HYPRE_ParVectorGetValues(vector,(int)indices_.size(),
            const_cast<HYPRE_BigInt*>(indices_.data()),result.data()),"Schur vector read");
        return result;
    }
    HYPRE_ParVector par() const{return par_;}
private:
    std::vector<HYPRE_BigInt> indices_;
    HYPRE_IJVector ij_=nullptr;
    HYPRE_ParVector par_=nullptr;
};

/// @brief 每次系数刷新时更新 Schur 响应，稀疏模式不变则复用 IJ 矩阵结构。
class SchurPreconditioner {
public:
    explicit SchurPreconditioner(FDM::LinearSolverConfig config):config_(config){}
    ~SchurPreconditioner(){clearPressure();}
    void update(const SparseSystem& system) {
        int ranks=0; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        if(ranks!=1)throw std::runtime_error("kktBlockSchur currently requires one MPI rank.");
        algebra_.assemble(system,config_.schurDenseLimit);
        const auto pressure=algebra_.pressureSystem();
        std::vector<std::vector<std::int64_t>> pattern;
        for(const auto& row:pressure.rows)pattern.push_back(row.columns);
        if(!matrix_ || pattern!=pattern_) {
            clearPressure(); pattern_=pattern;
            const int n=(int)pressure.rows.size();
            hypreCheck(HYPRE_IJMatrixCreate(MPI_COMM_WORLD,0,n-1,0,n-1,&matrix_),"Schur matrix create");
            hypreCheck(HYPRE_IJMatrixSetObjectType(matrix_,HYPRE_PARCSR),"Schur matrix type");
            b_=std::make_unique<SchurVector>(n); x_=std::make_unique<SchurVector>(n);
            hypreCheck(HYPRE_BoomerAMGCreate(&amg_),"Schur AMG create");
            configureAMGCycle(amg_,config_);
        }
        hypreCheck(HYPRE_IJMatrixInitialize(matrix_),"Schur matrix update");
        for(const auto& row:pressure.rows) {
            HYPRE_Int count=(int)row.columns.size();
            HYPRE_BigInt index=row.globalRow;
            std::vector<HYPRE_BigInt> columns(row.columns.begin(),row.columns.end());
            hypreCheck(HYPRE_IJMatrixSetValues(matrix_,1,&count,&index,columns.data(),row.values.data()),"Schur coefficients");
        }
        hypreCheck(HYPRE_IJMatrixAssemble(matrix_),"Schur matrix assemble");
        hypreCheck(HYPRE_IJMatrixGetObject(matrix_,reinterpret_cast<void**>(&parMatrix_)),"Schur matrix object");
        hypreCheck(HYPRE_BoomerAMGSetup(amg_,parMatrix_,b_->par(),x_->par()),"Schur AMG setup");
        algebra_.prepareConstraint([this](const auto& v){return pressureInverse(v);});
        if(!output_ || size_!=(int)system.rows.size()) {
            size_=(int)system.rows.size(); output_=std::make_unique<SchurVector>(size_);
        }
    }
    static HYPRE_Int setup(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector){return 0;}
    static HYPRE_Int apply(HYPRE_Solver context,HYPRE_ParCSRMatrix,HYPRE_ParVector rhs,HYPRE_ParVector solution) {
        auto& self=*reinterpret_cast<SchurPreconditioner*>(context);
        try {
            self.error_.clear();
            const auto result=self.algebra_.apply(self.output_->read(rhs),
                [&self](const auto& v){return self.pressureInverse(v);});
            self.output_->write(result);
            hypreCheck(HYPRE_ParVectorCopy(self.output_->par(),solution),"Schur output copy");
            return 0;
        } catch(const std::exception& error) {
            self.error_=error.what(); return HYPRE_ERROR_GENERIC;
        }
    }
    const std::string& error() const{return error_;}
private:
    std::vector<double> pressureInverse(const std::vector<double>& rhs) {
        b_->write(rhs);
        hypreCheck(HYPRE_ParVectorSetConstantValues(x_->par(),0),"Schur zero guess");
        hypreCheck(HYPRE_BoomerAMGSolve(amg_,parMatrix_,b_->par(),x_->par()),"Schur AMG cycle");
        return x_->read(x_->par());
    }
    void clearPressure() {
        if(amg_)HYPRE_BoomerAMGDestroy(amg_); amg_=nullptr;
        if(matrix_)HYPRE_IJMatrixDestroy(matrix_); matrix_=nullptr;
        parMatrix_=nullptr; b_.reset(); x_.reset();
    }
    FDM::LinearSolverConfig config_;
    BlockSchur algebra_;
    HYPRE_IJMatrix matrix_=nullptr;
    HYPRE_ParCSRMatrix parMatrix_=nullptr;
    HYPRE_Solver amg_=nullptr;
    std::unique_ptr<SchurVector> b_,x_,output_;
    std::vector<std::vector<std::int64_t>> pattern_;
    int size_=0;
    std::string error_;
};
} // namespace SF::LinearAlgebra
