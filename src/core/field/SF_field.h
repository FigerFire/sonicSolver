/// @file SF_field.h
/// @brief canonical 场数据、索引与状态存储实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/


#pragma once
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include "SF_valueTypes.h"
#include "SF_fieldComponents.h"
#include "SF_fluidStateModel.h"
#include "core/state/SF_thermodynamicStateCache.h"
#include <map>
#include <memory>

// memory layout 
// 
namespace SF {
    /// @brief 主守恒数组的连续存储次序。
    enum class FieldStorageLayout {
        AoS, ///< cell-major: Q[cell][variable]
        SoA  ///< variable-major: Q[variable][cell]
    };

    class Field {
    private:
        // 存储 RHO, RU, RV, RW, E
        std::vector<double> data; 
        
        // 维度信息 (包含虚胞)
        int  mx, my, mz, ng;
        int nVar_ = 5;
        int total_size;
        FieldStorageLayout storageLayout_ = FieldStorageLayout::AoS;
        std::shared_ptr<const Physics::FluidStateModel::Model> stateModel_;

        State::ThermodynamicStateCache thermodynamicCache_;

        GeometryStorage geometry_;
        BoundaryMetadata boundary_;
        /// @brief deferred IBM 分类标记（与 CellFlag 同类）；IBM 几何本体已迁出。
        std::vector<unsigned char> ibmFluidMask_;

        // 计算一维数组索引的辅助函数，必须加 const
        // 存储顺序是 [k][j][i][v]
        int calculateIndex(int i, int j, int k, int v) const {
             return static_cast<int>(idxV(i, j, k, v));
        }
           
        /// --- 私有索引辅助函数 ---
        inline size_t idxV(int i, int j, int k, int v) const {
            const size_t cell = idxS(i, j, k);
            return storageLayout_ == FieldStorageLayout::SoA
                ? (size_t)v * (size_t)total_size + cell
                : cell * (size_t)nVar_ + (size_t)v;
        }

        inline size_t idxS(int i, int j, int k) const {
            return (size_t)(k * my + j ) * mx + i;
        }

        /// @brief canonical face geometry 的平铺索引；不用于 solver workspace。
        inline size_t idxFace(int dir, int i, int j, int k) const {
            return (size_t)dir * (size_t)total_size + idxS(i, j, k);
        }

    public:
        /// @brief 默认构造,不分配内存。
        explicit Field(
            FieldStorageLayout layout = FieldStorageLayout::AoS);

        /// @brief 在分配前选择 AoS 或 SoA；已分配后禁止隐式重排。
        void setStorageLayout(FieldStorageLayout layout);
        FieldStorageLayout storageLayout() const { return storageLayout_; }

        /// @brief 返回连续主状态底层视图；布局由 storageLayout() 描述。
        double* conservativeData() { return data.data(); }
        const double* conservativeData() const { return data.data(); }

        /// @brief SoA 模式下返回一个变量的连续数组；AoS 模式调用会失败。
        double* componentData(int variable);
        const double* componentData(int variable) const;

        /// @brief 初始化Field的维度和存储空间。
        /// @param nx, ny, nz 内部物理网格点数(不含虚胞)。
        /// @param ghost 虚胞层数(每侧)。
        void setup(int nx, int ny, int nz, int ghost, int nVar = 0);

        /// @brief 绑定控制守恒变量布局和热力学闭合的 FluidStateModel。
        /// 已分配 Field 的变量数必须与方程组匹配；需要切换布局时先调用
        /// resizeConservedVariables，再写入新的守恒状态。
        void setStateModel(std::shared_ptr<const Physics::FluidStateModel::Model> equations);
        const std::shared_ptr<const Physics::FluidStateModel::Model>& stateModel() const { return stateModel_; }
        bool hasStateModel() const { return static_cast<bool>(stateModel_); }

        /// @brief 保留网格、集合和分类标记，重分配可变守恒变量存储。
        /// 旧守恒量不会被隐式投影或复制；调用者必须立即按新 FluidStateModel 重建状态。
        void resizeConservedVariables(int nVar);

        /// @brief 读取一个点的 EOS-aware 派生热力学状态。
        Physics::FluidStateModel::ThermodynamicState thermodynamicState(int i, int j, int k) const;

        /// @brief 声明守恒时间层已写入，使 EOS 派生缓存失效。
        /// @details 时间积分、边界、IBM、halo 和压力修正完成写入后必须调用。
        void invalidateThermodynamicCache();

        /// @brief 当前守恒状态版本，用于诊断缓存生命周期。
        std::uint64_t stateVersion() const { return thermodynamicCache_.version(); }


        /// @brief 从一维平铺索引逆向计算(i,j,k)。
        /// @param idx 0-based一维平铺索引。
        /// @param i, j, k 输出的三维索引。
        void getIJK(int idx, int& i, int& j, int& k) const;

        // ── 集合索引访问 ──

        /// @brief 设置边界集合点映射(整体替换)。
        /// @param sets key=集合名, value=Field内部索引列表。
        void setBoundarySets(const std::map<std::string, std::vector<int>>& sets);

        /// @brief 获取指定名称的边界集合点列表。
        /// @param name 集合名(如"Left", "Right", "Wall")。
        /// @return 索引列表的const引用; 不存在时返回空vector。
        const std::vector<int>& getSet(const std::string& name) const;

        // ── 维度访问器 ──

        /// @brief 含虚胞的总X方向网格点数。
        int MX() const { return mx; }
        /// @brief 含虚胞的总Y方向网格点数。
        int MY() const { return my; }
        /// @brief 含虚胞的总Z方向网格点数。
        int MZ() const { return mz; }
        /// @brief 含虚胞的总网格点数 (=MX*MY*MZ)。
        int TotalSize() const { return total_size; }
        /// @brief 虚胞层数。
        int NG() const { return ng; }
        /// @brief 当前 FluidStateModel 的守恒变量数量。
        int NVar() const { return nVar_; }

        /// @brief 纯内部网格X方向点数(不含虚胞)。
        int NX() const { return mx - 2 * ng; }
        /// @brief 纯内部网格Y方向点数(不含虚胞)。
        int NY() const { return my - 2 * ng; }
        /// @brief 纯内部网格Z方向点数(不含虚胞)。
        int NZ() const { return mz - 2 * ng; }

        /// @brief 清空由BC set生成的求解边界mask。
        void clearSolverBoundaryMask();

        /// @brief 将一个BC集合标记为真实物理边界求解边界。
        /// @param name BC配置中引用的集合名。
        void markSolverBoundarySet(const std::string& name);

        /// @brief 判断平铺索引是否属于指定集合。
        /// @param name 集合名。
        /// @param idx 0-based平铺索引。
        /// @return idx存在于集合中时返回true。
        bool isInSet(const std::string& name, int idx) const;

        /// @brief 将某个坐标方向两端物理面标记为真实物理边界。
        /// @param axis 0=XI, 1=ETA, 2=ZETA。
        void markSolverBoundaryAxis(int axis);

        /// @brief 判断给定点是否属于真实物理边界求解边界。
        /// @param i i方向全局存储索引。
        /// @param j j方向全局存储索引。
        /// @param k k方向全局存储索引。
        /// @return 属于BC set声明的真实物理边界时返回true。
        bool isSolverBoundaryPoint(int i, int j, int k) const;

        // ── 索引计算 ──

        /// @brief 计算(i,j,k)在平铺数组中的1D索引。
        /// @param i, j, k 三维全局索引(含虚胞)。
        /// @return 0-based平铺索引。
        inline int getIdx(int i, int j, int k) const { return (k * my + j) * mx + i; }

        /// @brief 获取可写的边界集合索引列表(不存在时会创建)。
        /// @param name 集合名。
        /// @return 索引列表的可写引用。
        std::vector<int>& getSetWritable(const std::string& name) {
            return boundary_.sets[name];
        }

        /// @brief 获取所有边界集合的映射(可写)。
        /// @return boundarySets的可写引用。
        std::map<std::string, std::vector<int>>& getAllSets() {
            return boundary_.sets;
        }
        /// @brief 获取所有边界集合的映射(只读)。
        /// @return boundarySets的const引用。
        const std::map<std::string, std::vector<int>>& getAllSets() const {
            return boundary_.sets;
        }

        // ── 单元标记 ──

        /// @brief 获取/设置(i,j,k)处网格单元的类型标记。
        /// @param i, j, k 三维全局索引。
        /// @return 单元类型(FLUID_CELL, IBM_GHOST_CELL, SOLID_CELL)的引用。
        inline int& CellFlag(int i, int j, int k) {
            return boundary_.cellType[idxS(i, j, k)];
        }
        /// @brief 获取(i,j,k)处网格单元的类型标记(只读)。
        inline int CellFlag(int i, int j, int k) const {
            return boundary_.cellType[idxS(i, j, k)];
        }

        /// @brief 全网格单元类型重置为指定标记。
        /// @param marker 目标标记, 默认FLUID_CELL。
        void clearCellFlags(int marker = FLUID_CELL);

        /// @brief 清空分区/多块接口 halo 标记；Field 不执行通信。
        void clearCommunicationHaloMask();

        /// @brief 设置一个存储点是否为已映射的分区/多块接口 halo。
        inline void setCommunicationHalo(int i, int j, int k, bool enabled) {
            boundary_.communicationHaloMask[idxS(i, j, k)] =
                enabled ? 1 : 0;
        }

        /// @brief 查询存储点是否为已映射的分区/多块接口 halo。
        inline bool isCommunicationHalo(int i, int j, int k) const {
            return boundary_.communicationHaloMask[idxS(i, j, k)] != 0;
        }

        /// @brief 设置一个 Eulerian 点的 GlobalDof 及唯一 owner。
        inline void setGlobalDofOwnership(
                int i, int j, int k, int globalDofId,
                int ownerRank, bool isOwner) {
            const size_t id = idxS(i, j, k);
            boundary_.globalDofId[id] = globalDofId;
            boundary_.globalDofOwnerRank[id] = ownerRank;
            boundary_.globalDofOwnerMask[id] = isOwner ? 1 : 0;
        }

        /// @brief 返回 Eulerian 点的 GlobalDof；未注册的 ghost 返回 -1。
        inline int globalDofId(int i, int j, int k) const {
            return boundary_.globalDofId[idxS(i, j, k)];
        }

        /// @brief 返回 GlobalDof 的唯一 MPI owner rank。
        inline int globalDofOwnerRank(int i, int j, int k) const {
            return boundary_.globalDofOwnerRank[idxS(i, j, k)];
        }

        /// @brief 当前结构 patch 副本是否为 GlobalDof canonical owner。
        inline bool isGlobalDofOwner(int i, int j, int k) const {
            return boundary_.globalDofOwnerMask[idxS(i, j, k)] != 0;
        }

        inline bool IBMFluidMask(int i, int j, int k) {
            return ibmFluidMask_[idxS(i, j, k)] != 0;
        }
        inline bool IBMFluidMask(int i, int j, int k) const {
            return ibmFluidMask_[idxS(i, j, k)] != 0;
        }
        inline void setIBMFluidMask(int i, int j, int k, bool outsideFluid) {
            ibmFluidMask_[idxS(i, j, k)] = outsideFluid ? 1 : 0;
        }

        /// --- 1. 守恒变量访问 (data 数组) ---
        // Int 版本
        double& operator()(int i, int j, int k, int v) { return data[idxV(i, j, k, v)]; }
        double  operator()(int i, int j, int k, int v) const { return data[idxV(i, j, k, v)]; }

        /// 枚举版本 (VarIdx v: RHO, RU, RV, RW, E)
        double& operator()(int i, int j, int k, VarIdx v) { return (*this)(i, j, k, static_cast<int>(v)); }
        double  operator()(int i, int j, int k, VarIdx v) const { return (*this)(i, j, k, static_cast<int>(v)); }

        /// --- 2. 坐标访问 (x, y, z 数组) ---
        inline double& X(int i, int j, int k) { return geometry_.x[idxS(i, j, k)]; }
        inline double X(int i, int j, int k) const { return geometry_.x[idxS(i, j, k)]; }
        inline double& Y(int i, int j, int k) { return geometry_.y[idxS(i, j, k)]; }
        inline double Y(int i, int j, int k) const { return geometry_.y[idxS(i, j, k)]; }
        inline double& Z(int i, int j, int k) { return geometry_.z[idxS(i, j, k)]; }
        inline double Z(int i, int j, int k) const { return geometry_.z[idxS(i, j, k)]; }

        // --- 2.5 到最近壁面的有符号距离 ---
        inline double& wallDistance(int i, int j, int k) {
            return geometry_.wallDistance[idxS(i, j, k)];
        }
        inline double wallDistance(int i, int j, int k) const {
            return geometry_.wallDistance[idxS(i, j, k)];
        }
        inline void setWallDistance(int i, int j, int k, double d) {
            geometry_.wallDistance[idxS(i, j, k)] = d;
        }

        // --- 3. 雅可比访问 (jac 数组) ---
        inline double& Jac(int i, int j, int k) { return geometry_.jac[idxS(i, j, k)]; }
        inline double Jac(int i, int j, int k) const { return geometry_.jac[idxS(i, j, k)]; }

        // --- 4. 度规访问 (xi, et, ze 数组) ---
        // Xi 方向 (对应 i)
        inline double& XiX(int i, int j, int k) { return geometry_.xiX[idxS(i, j, k)]; }
        inline double XiX(int i, int j, int k) const { return geometry_.xiX[idxS(i, j, k)]; }
        inline double& XiY(int i, int j, int k) { return geometry_.xiY[idxS(i, j, k)]; }
        inline double XiY(int i, int j, int k) const { return geometry_.xiY[idxS(i, j, k)]; }
        inline double& XiZ(int i, int j, int k) { return geometry_.xiZ[idxS(i, j, k)]; }
        inline double XiZ(int i, int j, int k) const { return geometry_.xiZ[idxS(i, j, k)]; }

        // Eta 方向 (对应 j)
        inline double& EtX(int i, int j, int k) { return geometry_.etaX[idxS(i, j, k)]; }
        inline double EtX(int i, int j, int k) const { return geometry_.etaX[idxS(i, j, k)]; }
        inline double& EtY(int i, int j, int k) { return geometry_.etaY[idxS(i, j, k)]; }
        inline double EtY(int i, int j, int k) const { return geometry_.etaY[idxS(i, j, k)]; }
        inline double& EtZ(int i, int j, int k) { return geometry_.etaZ[idxS(i, j, k)]; }
        inline double EtZ(int i, int j, int k) const { return geometry_.etaZ[idxS(i, j, k)]; }

        // Zeta 方向 (对应 k)
        inline double& ZeX(int i, int j, int k) { return geometry_.zetaX[idxS(i, j, k)]; }
        inline double ZeX(int i, int j, int k) const { return geometry_.zetaX[idxS(i, j, k)]; }
        inline double& ZeY(int i, int j, int k) { return geometry_.zetaY[idxS(i, j, k)]; }
        inline double ZeY(int i, int j, int k) const { return geometry_.zetaY[idxS(i, j, k)]; }
        inline double& ZeZ(int i, int j, int k) { return geometry_.zetaZ[idxS(i, j, k)]; }
        inline double ZeZ(int i, int j, int k) const { return geometry_.zetaZ[idxS(i, j, k)]; }

        /// @brief 为接口半面写入唯一 canonical cofactor/Jacobian。
        void setCanonicalFaceMetrics(int dir, int i, int j, int k,
                                     const std::array<double, 4>& metrics);
        /// @brief 查询接口半面是否有 canonical geometry。
        bool hasCanonicalFaceMetrics(int dir, int i, int j, int k) const;
        /// @brief 读取接口半面的 canonical geometry。
        void canonicalFaceMetrics(int dir, int i, int j, int k,
                                  double metrics[4]) const;
        /// @brief 标记当前结构半面是否为 GlobalFace 的唯一计算 owner。
        void setCanonicalFaceOwner(
            int dir, int i, int j, int k, bool isOwner);
        /// @brief 查询当前结构半面是否负责计算 canonical 数值通量。
        bool isCanonicalFaceOwner(int dir, int i, int j, int k) const;


        // 在 Field 类中
        /// @brief 读取一个变量组：Vector3 由 v 指定第一个分量。
        /// @param v 变量组的第一个连续分量索引（例如 rhoU 组为 1）。
        /// @details Vector3 的布局由调用方通过 StateModel/StorageBinding
        ///          给出的分量起点决定；这里不再假定 RU/RV/RW 固定为 1/2/3，
        ///          否则任何非 1/2/3 的 vector 组都会被静默读写错位。
        template<typename T>
        T getVal(int i, int j, int k, int v) {
            if constexpr (std::is_same_v<T, Vector3>)
                return Vector3((*this)(i,j,k,v), (*this)(i,j,k,v+1), (*this)(i,j,k,v+2));
            else
                return (*this)(i,j,k,v);
        }

        template<typename T>
        void setVal(int i, int j, int k, int v, const T& val) {
            if constexpr (std::is_same_v<T, Vector3>) {
                (*this)(i,j,k,v)   = val.x;
                (*this)(i,j,k,v+1) = val.y;
                (*this)(i,j,k,v+2) = val.z;
            } else {
                (*this)(i,j,k,v) = val;
            }
        }
    };
}
