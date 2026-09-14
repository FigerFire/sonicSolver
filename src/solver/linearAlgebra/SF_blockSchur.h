#pragma once
/// @file SF_blockSchur.h
/// @brief 串行 KKT 块消元：速度对角近似、刚体稠密块及压力/乘子交叉耦合。
#include "SF_linearAlgebra.h"
#include "methods/math/discrete/SF_denseLU.h"
#include <functional>
#include <map>

namespace SF::LinearAlgebra {
/// @brief 后端无关的块三角分解；pressureInverse 是固定线性 AMG 周期。
class BlockSchur {
public:
    using Row=std::map<int,double>;
    using Inverse=std::function<std::vector<double>(const std::vector<double>&)>;
    void assemble(const SparseSystem& system,int denseLimit) {
        const int n=(int)system.rows.size();
        if(system.firstRow!=0 || system.globalSize!=n
           ||system.blockRoles.size()!=system.rows.size())
            throw std::runtime_error("kktBlockSchur requires serial, labelled KKT rows.");
        roles_=system.blockRoles; p_.clear(); l_.clear(); a_.clear(); solid_.clear();
        for(int i=0;i<n;++i) {
            if(roles_[i]==1)p_.push_back(i);
            else if(roles_[i]==3)l_.push_back(i);
            else if(roles_[i]==2 || roles_[i]==4) {
                if(roles_[i]==4)solid_.push_back((int)a_.size());
                a_.push_back(i);
            } else throw std::runtime_error("kktBlockSchur unknown row role.");
        }
        if(p_.empty() || l_.empty() || (int)l_.size()>denseLimit)
            throw std::runtime_error("kktBlockSchur missing blocks or schurDenseLimit exceeded.");
        reduced_=p_; reduced_.insert(reduced_.end(),l_.begin(),l_.end());
        std::vector<int> ai(n,-1),ri(n,-1);
        for(int i=0;i<(int)a_.size();++i)ai[a_[i]]=i;
        for(int i=0;i<(int)reduced_.size();++i)ri[reduced_[i]]=i;
        diagonal_.assign(a_.size(),0); b_.assign(a_.size(),{});
        c_.assign(reduced_.size(),{}); schur_.assign(reduced_.size(),{});
        const int ns=(int)solid_.size();
        std::vector<int> si(a_.size(),-1);
        for(int s=0;s<ns;++s)si[solid_[s]]=s;
        std::vector<double> solidMatrix(ns*ns,0);
        for(int row=0;row<n;++row) {
            const auto& input=system.rows[row];
            for(size_t j=0;j<input.columns.size();++j) {
                const int col=(int)input.columns[j]; const double v=input.values[j];
                if(ai[row]>=0) {
                    if(ri[col]>=0)b_[ai[row]][ri[col]]+=v;
                    else if(row==col)diagonal_[ai[row]]+=v;
                    if(roles_[row]==4 && roles_[col]==4)
                        solidMatrix[si[ai[row]]*ns+si[ai[col]]]+=v;
                    // 速度非对角项只在预条件器里采用对角近似；完整外层矩阵保留。
                    if(ai[col]>=0 && roles_[row]!=roles_[col] && v!=0)
                        throw std::runtime_error("kktBlockSchur unexpected primal cross block.");
                } else if(ai[col]>=0)c_[ri[row]][ai[col]]+=v;
                else schur_[ri[row]][ri[col]]+=v;
            }
        }
        for(size_t i=0;i<a_.size();++i)
            if(roles_[a_[i]]==2 && (!std::isfinite(diagonal_[i]) || diagonal_[i]<=0))
                throw std::runtime_error("kktBlockSchur non-positive velocity diagonal.");
        if(ns)solidLU_.factor(solidMatrix,ns);
        inverseB_=b_;
        for(size_t i=0;i<a_.size();++i)if(roles_[a_[i]]==2)
            for(auto& item:inverseB_[i])item.second/=diagonal_[i];
        std::map<int,bool> solidColumns;
        for(int s:solid_)for(auto item:b_[s])solidColumns[item.first]=true;
        for(auto column:solidColumns) {
            std::vector<double> rhs(ns,0);
            for(int s=0;s<ns;++s)rhs[s]=b_[solid_[s]][column.first];
            const auto response=solidLU_.solve(rhs);
            for(int s=0;s<ns;++s)inverseB_[solid_[s]][column.first]=response[s];
        }
        for(size_t row=0;row<c_.size();++row)
            for(auto entry:c_[row])for(auto col:inverseB_[entry.first])
                schur_[row][col.first]-=entry.second*col.second;
    }
    SparseSystem pressureSystem() const {
        SparseSystem result; const int np=(int)p_.size();
        result.globalSize=np; result.lastRow=np-1;
        result.rhs.assign(np,0);
        for(int i=0;i<np;++i) {
            SparseRow row; row.globalRow=i;
            for(auto entry:schur_[i])if(entry.first<np && entry.second!=0) {
                row.columns.push_back(entry.first); row.values.push_back(entry.second);
            }
            result.rows.push_back(std::move(row));
        }
        return result;
    }
    void prepareConstraint(const Inverse& pressureInverse) {
        const int np=(int)p_.size(), nl=(int)l_.size();
        pressureResponse_.assign(nl,{});
        std::vector<double> dense(nl*nl,0);
        for(int col=0;col<nl;++col) {
            std::vector<double> rhs(np,0);
            for(int i=0;i<np;++i)rhs[i]=value(schur_[i],np+col);
            pressureResponse_[col]=pressureInverse(rhs);
            for(int row=0;row<nl;++row) {
                double v=value(schur_[np+row],np+col);
                for(auto e:schur_[np+row])if(e.first<np)
                    v-=e.second*pressureResponse_[col][e.first];
                dense[row*nl+col]=v;
            }
        }
        constraintLU_.factor(dense,nl);
    }
    std::vector<double> apply(const std::vector<double>& rhs,const Inverse& pressureInverse) const {
        const int np=(int)p_.size(),nl=(int)l_.size();
        if(rhs.size()!=roles_.size())throw std::runtime_error("BlockSchur RHS size mismatch.");
        std::vector<double> xA(a_.size(),0), r(reduced_.size());
        for(size_t i=0;i<a_.size();++i)if(roles_[a_[i]]==2)xA[i]=rhs[a_[i]]/diagonal_[i];
        if(!solid_.empty()) {
            std::vector<double> s;
            for(int i:solid_)s.push_back(rhs[a_[i]]);
            s=solidLU_.solve(s);
            for(size_t i=0;i<solid_.size();++i)xA[solid_[i]]=s[i];
        }
        for(size_t i=0;i<r.size();++i) {
            r[i]=rhs[reduced_[i]];
            for(auto e:c_[i])r[i]-=e.second*xA[e.first];
        }
        auto pressure=pressureInverse(std::vector<double>(r.begin(),r.begin()+np));
        std::vector<double> lambda(nl);
        for(int i=0;i<nl;++i) {
            lambda[i]=r[np+i];
            for(auto e:schur_[np+i])if(e.first<np)lambda[i]-=e.second*pressure[e.first];
        }
        lambda=constraintLU_.solve(lambda);
        for(int i=0;i<np;++i)for(int j=0;j<nl;++j)
            pressure[i]-=pressureResponse_[j][i]*lambda[j];
        std::vector<double> y=pressure; y.insert(y.end(),lambda.begin(),lambda.end());
        std::vector<double> result(rhs.size(),0);
        for(size_t i=0;i<a_.size();++i) {
            for(auto e:inverseB_[i])xA[i]-=e.second*y[e.first];
            result[a_[i]]=xA[i];
        }
        for(size_t i=0;i<y.size();++i)result[reduced_[i]]=y[i];
        return result;
    }
private:
    static double value(const Row& row,int col) {
        auto it=row.find(col); return it==row.end()?0:it->second;
    }
    std::vector<int> roles_,p_,l_,a_,solid_,reduced_;
    std::vector<double> diagonal_;
    std::vector<Row> b_,c_,inverseB_,schur_;
    std::vector<std::vector<double>> pressureResponse_;
    Math::DenseLU solidLU_,constraintLU_;
};
} // namespace SF::LinearAlgebra
