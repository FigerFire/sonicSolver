/// @file SF_phaseProperties.h
/// @brief 多相配置、物性与模型一致性校验实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_valueTypes.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set 符号约定中的相角色。
///
/// 当前约定: liquid 为 phi > 0, gas 为 phi < 0。
enum class PhaseRole {
    Unknown,
    Liquid,
    Gas
};

/// @brief 相内组分；Y 是该相内质量分数，而非流域总体分数。
struct SpeciesProperties {
    std::string name;
    double massFraction = -1.0; ///< 均质 FluidStateModel 初值 Y；必须显式且同相求和为一。
};

/// @brief 固定压力下的温度物性定律。
///
/// `polynomialTemperature` 使用
/// `f(T)=sum(coefficients[i]*(T-referenceTemperature)^i)`。多项式必须声明
/// 有效温区，求解状态越界时 fail-fast，不能静默外推。
struct TemperaturePropertyLaw {
    std::string model = "constant";
    double referenceTemperature = 0.0;
    double minimumTemperature = 0.0;
    double maximumTemperature = 0.0;
    std::vector<double> coefficients;
};

/// @brief 单个相的物性和符号角色。
struct PhaseProperties {
    std::string name;              ///< 相名, 如 water/air。
    PhaseRole role = PhaseRole::Unknown;
    bool roleDeclared = false;     ///< role/明确 EOS 是否由用户声明；Level Set 不按相名猜测。
    std::string thermoModel;       ///< perfectLiquid / perfectGas。
    double density = 0.0;          ///< rho0/rho；未配置时为0, 不参与隐式兜底。
    double viscosity = 0.0;        ///< 动力黏度；未配置时为0, 不参与隐式兜底。
    double specificHeat = 0.0;     ///< 定压比热 Cp；压力基相焓和热输运使用。
    double Cv = 0.0;               ///< 定容比热；密度基 rhoE 和 EOS 使用。
    double e0 = 0.0;               ///< T=T_ref 时的参考内能 (J/kg)。
    double T_ref = 273.15;         ///< 参考温度 (K)；默认 273.15 K。
    double gamma = 0.0;            ///< EOS 比热比；不使用全局固定 gamma。
    double gasConstant = 0.0;      ///< PerfectGas EOS 气体常数。
    double pInfinity = 0.0;        ///< StiffenedGas EOS 刚化压力。
    double volumeFraction = -1.0;  ///< 均质 FluidStateModel 初值相体积分数；<0 表示未声明。
    std::vector<SpeciesProperties> species; ///< 空列表退化为一个名为 bulk 的组分。
    double thermalConductivity = 0.0; ///< 导热系数；需要热通量闭式关系时显式给出。
    TemperaturePropertyLaw densityLaw; ///< rho(T)；未声明时使用 EOS/常密度。
    TemperaturePropertyLaw viscosityLaw; ///< 动力黏度 mu(T)。
    TemperaturePropertyLaw specificHeatLaw; ///< Cp(T)，同时定义显焓积分。
    TemperaturePropertyLaw thermalConductivityLaw; ///< k(T)。
    std::string diameterModel = "none"; ///< none/constant。
    double diameter = 0.0;             ///< 分散相代表直径 (m)。
    double residualAlpha = 0.0;        ///< 闭式模型的显式残余体积分数。
    bool initialAlphaDeclared = false;
    bool initialDensityDeclared = false;
    bool initialVelocityDeclared = false;
    bool initialTemperatureDeclared = false;
    Vector3 initialVelocity;
    double initialTemperature = std::numeric_limits<double>::quiet_NaN();
    std::vector<BCSetting<double>> alphaInitial;
    std::vector<BCSetting<double>> alphaBoundary;
    std::vector<BCSetting<double>> densityInitial;
    std::vector<BCSetting<double>> densityBoundary;
    std::vector<BCSetting<Vector3>> velocityInitial;
    std::vector<BCSetting<Vector3>> velocityBoundary;
    std::vector<BCSetting<double>> temperatureInitial;
    std::vector<BCSetting<double>> temperatureBoundary;
};

/// @brief 一对 Eulerian 相之间的闭式配置。
///
/// 每个 pair 独立声明连续相、分散相和闭式模型；PhaseSystem 不再保存唯一的
/// continuous/dispersed 角色，因此三相及更多相时可以注册任意无重复相对。
struct PhasePairModelOptions {
    std::string name;
    std::string continuousPhase;
    std::string dispersedPhase;
    std::string dragModel = "SchillerNaumann";
    std::string liftModel = "constantCoefficient";
    std::string virtualMassModel = "constantCoefficient";
    std::string heatTransferModel = "RanzMarshall";
    std::string wallLubricationModel = "none";
    std::string wallPatch;
    std::string turbulentDispersionModel = "none";
    std::string surfaceTensionModel = "none";
    double surfaceTension = 0.0;
    double particleDiameter = 0.0;
    double liftCoefficient = 0.0;
    double virtualMassCoefficient = 0.5;
    double wallLubricationC1 = -0.01;
    double wallLubricationC2 = 0.05;
    double turbulentDispersionCoefficient = 0.0;
    double turbulentKinematicViscosity = 0.0;
    double turbulentSchmidtNumber = 0.0;
};

/// @brief Eulerian-Eulerian PhaseSystem 的配置。
struct EulerianEulerianOptions {
    std::vector<std::string> phaseNames;
    std::string referencePhase;
    bool axisymmetric = false; ///< 2D r-z 控制体采用 2*pi*r 几何权重。
    int radialCoordinate = 0;  ///< 0/1/2 对应 x/y/z。
    std::vector<PhasePairModelOptions> phasePairs;
};

/// @brief 多相求解器家族。
enum class MultiPhaseSolverKind {
    Auto,
    DensityBased,
    PressureBased
};

/// @brief mixture 模型输入结构。
struct MixtureOptions {
    std::vector<std::string> phaseNames; ///< phases (...) 中的相名。
    std::vector<std::string> phaseFractionFields; ///< alpha.liquid 等相分数字段。
};

/// @brief density-based 扩展守恒量布局说明。
struct ConservativeLayout {
    std::vector<std::string> variables; ///< rho/rhoU/rhoE/rhoYi/rhoAlpha 顺序。
};

/// @brief 网格 set 到相名的初始化映射。
struct SetPhaseAssignment {
    std::string setName;           ///< Field 中已有的 set 名。
    std::string phaseName;         ///< 目标相名。
};

/// @brief Level Set 标量场的 IC/BC 条件。
struct LevelSetScalarCondition {
    std::string setName;           ///< Field 中已有的 set 名。
    BCType type = FIXED_VALUE;     ///< FIXED_VALUE/ZERO_GRADIENT/EMPTY 等。
    bool typeDeclared = false;     ///< 条件类型是否由外部输入显式给出。
    double value = std::numeric_limits<double>::quiet_NaN(); ///< phi 数值；phi>0 为液相。
};

/// @brief Level Set 平面 signed-distance 初始条件。
struct LevelSetPlaneInitializer {
    std::string setName;             ///< 需要写入 signed-distance 的 set。
    Vector3 point = Vector3(
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()); ///< 平面上一点。
    Vector3 normal = Vector3(
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()); ///< 指向 phi>0 一侧的法向。
    double scale = std::numeric_limits<double>::quiet_NaN(); ///< signed-distance 乘子。
};

/// @brief Level Set 球面/二维圆形 signed-distance 初始条件。
struct LevelSetSphereInitializer {
    std::string setName;             ///< 需要写入 signed-distance 的 set。
    Vector3 center = Vector3(
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()); ///< 球心或二维圆心。
    double radius = std::numeric_limits<double>::quiet_NaN(); ///< 半径（m）。
    double scale = std::numeric_limits<double>::quiet_NaN(); ///< signed-distance 乘子。
    std::string insidePhase;         ///< 球内相名；用于显式确定 phi 符号。
};

/// @brief Level Set 子模块配置。
struct LevelSetOptions {
    int advectionOrder = 0;         ///< 必填：对流 HJ-WENO 阶数, 3/5/7。
    double wenoEpsilon = std::numeric_limits<double>::quiet_NaN(); ///< 必填：非线性权重正则量。
    double wenoPower = std::numeric_limits<double>::quiet_NaN(); ///< 必填：非线性权重指数。
    bool advectFluidCellsOnly = false; ///< 是否只推进 FLUID_CELL。
    bool advectFluidCellsOnlyDeclared = false; ///< 上述开关是否显式输入。
    int reinitializationOrder = 0;  ///< 必填：重初始化 HJ-WENO 阶数, 3/5/7。
    int reinitializationSteps = -1; ///< 必填：重初始化伪时间步数；0 表示关闭。
    double pseudoTimeStep = std::numeric_limits<double>::quiet_NaN(); ///< 必填：伪时间步长；关闭时显式写0。
    double signSmoothingFactor = std::numeric_limits<double>::quiet_NaN(); ///< 必填：S(phi0) 平滑宽度相对最小网格尺度的倍率。
    double geometryGradientTolerance = std::numeric_limits<double>::quiet_NaN(); ///< 必填：法向计算最小梯度。
    double interfaceThickness = std::numeric_limits<double>::quiet_NaN(); ///< 必填：Heaviside/delta 半厚度；0 为锐界面物性。
    double narrowBandWidth = std::numeric_limits<double>::quiet_NaN(); ///< 必填：窄带半宽；0 表示关闭。
    std::string surfaceTensionModel; ///< 必填：none/CSF/ghostFluid。
    bool surfaceTensionModelDeclared = false; ///< 上述模型是否由用户显式输入。
    double surfaceTension = std::numeric_limits<double>::quiet_NaN(); ///< 必填：表面张力系数（N/m）。
    bool massCorrection = false;    ///< 质量修正开关；当前显式启用会 fail-fast。
    bool massCorrectionDeclared = false; ///< 上述开关是否显式输入。
};

/// @brief mixture-like 单速度模型中的 alpha 标量配置。
struct AlphaOptions {
    std::string fieldName = "alpha"; ///< VTK/诊断中使用的 alpha 名称。
    std::string phaseName = "liquid"; ///< alpha=1 对应的相名。
    double defaultValue = 0.0;        ///< 未被 IC 覆盖区域的初始 alpha。
    bool transportEnabled = true;     ///< 是否随单速度场输运 alpha。
    bool diffusionEnabled = false;    ///< 是否启用常扩散项。
    double diffusivity = 0.0;         ///< 常扩散系数。
    bool bounded = true;              ///< 是否启用 [lower, upper] 诊断/投影。
    double lowerBound = 0.0;          ///< alpha 下界。
    double upperBound = 1.0;          ///< alpha 上界。
    std::string boundMode = "diagnose"; ///< none/diagnose/project(clip)。
    std::vector<BCSetting<double>> initialConditions; ///< 0/alpha 初值。
    std::vector<BCSetting<double>> boundaryConditions;///< 0/alpha 边界。
};

/// @brief 单温度热输运标量配置。
///
/// 当前热模型是独立的温度标量方程, 不把 `T` 从 Euler 总能量的理想气体闭合中
/// 反推出来。相变潜热源当前接入 Euler 总能量源项; 后续若要严格做混合物焓
/// 守恒耦合, 仍需要新增 `h` 或 `rhoH` 的保守输运路径, 并统一刷新温度闭合。
struct TemperatureOptions {
    bool enabled = false;             ///< 是否启用温度标量方程。
    std::string fieldName = "T";      ///< VTK/诊断中使用的温度名。
    double defaultValue = std::numeric_limits<double>::quiet_NaN(); ///< 显式温度默认值。
    bool transportEnabled = true;     ///< 是否随单速度场输运温度。
    bool diffusionEnabled = false;    ///< 是否启用常扩散项。
    double diffusivity = 0.0;         ///< 温度常扩散系数。
    bool bounded = false;             ///< 是否启用温度上下界诊断/投影。
    double lowerBound = 0.0;          ///< 温度下界。
    double upperBound = 1.0e300;      ///< 温度上界。
    std::string boundMode = "diagnose"; ///< none/diagnose/project(clip)。
    std::vector<BCSetting<double>> initialConditions; ///< 0/T 初值。
    std::vector<BCSetting<double>> boundaryConditions;///< 0/T 标量边界。
};

/// @brief 简单过热相变源项配置。
struct PhaseChangeOptions {
    bool enabled = false;             ///< 是否启用相变源项。
    std::vector<std::string> phaseNames; ///< 当前相变涉及的有序相对。
    std::string model = "none";       ///< none / superheat。
    double saturationTemperature = 0.0; ///< 饱和温度阈值。
    double saturationPressure = 0.0;  ///< 饱和压力；空化/HKS模型需要。
    double evaporationCoefficient = 0.0; ///< T>Tsat 时的 alpha 源系数。
    double condensationCoefficient = 0.0; ///< T<Tsat 时的 alpha 源系数。
    double latentHeat = 0.0;          ///< 正潜热；蒸发吸热、冷凝放热。
    bool energyCoupling = true;       ///< 相变必须耦合能量，false 会 fail-fast。
    std::map<std::string, double> coefficients; ///< 模型专用显式系数。
    std::map<std::string, std::string> selections; ///< 模型专用显式选项。
    std::vector<WallHeatSetting> wallBoiling; ///< RPI-like 壁面沸腾热流设置。
};

/// @brief 多相守恒诊断量。
///
/// 这些量只做诊断，不反向修正解。`phaseIndicatorIntegral` 是
/// `alpha` 或 Level Set Heaviside 在物理体积上的积分；`phaseMass`
/// 进一步乘以被跟踪相密度。误差相对初始化基准计算。
struct ConservationDiagnostics {
    bool initialized = false;          ///< 是否已有初始化基准。
    double referencePhaseIndicator = 0.0;
    double currentPhaseIndicator = 0.0;
    double phaseIndicatorError = 0.0;
    double phaseIndicatorRelativeError = 0.0;
    double referencePhaseMass = 0.0;
    double currentPhaseMass = 0.0;
    double phaseMassError = 0.0;
    double phaseMassRelativeError = 0.0;
    double referenceFlowMass = 0.0;
    double currentFlowMass = 0.0;
    double flowMassError = 0.0;
    double flowMassRelativeError = 0.0;
};

/// @brief 多相/热物理模型配置对象, 由 IO 层解析后传给 physics 层。
struct MultiPhaseConfig {
    bool enabled = false;
    std::string type = "none";                 ///< level set / mixture / thermal / none。
    MultiPhaseSolverKind solver = MultiPhaseSolverKind::Auto; ///< auto/densityBased/pressureBased。
    std::string defaultPhase;                   ///< Level Set 必填：未被 set 覆盖的相。
    double defaultSignedDistance = std::numeric_limits<double>::quiet_NaN(); ///< Level Set 必填：set 初始化符号场幅值。
    double initialPressure = 0.0;               ///< homogeneous EOS 初值共压；必须显式给出。
    MixtureOptions mixture;
    ConservativeLayout conservativeLayout;
    std::vector<PhaseProperties> phases;
    std::vector<SetPhaseAssignment> setPhases;  ///< phaseProperties 中的 set->phase 初始化。
    std::vector<LevelSetPlaneInitializer> phiPlaneInitializers;
    std::vector<LevelSetSphereInitializer> phiSphereInitializers;
    std::vector<LevelSetScalarCondition> phiInitialConditions;
    std::vector<LevelSetScalarCondition> phiBoundaryConditions;
    LevelSetOptions levelSet;
    AlphaOptions alpha;
    TemperatureOptions temperature;
    PhaseChangeOptions phaseChange;
    EulerianEulerianOptions eulerianEulerian;
};

/// @brief 规整多相模型类型名。
/// @param type 原始类型名。
/// @return 去掉空格/下划线/连字符并转小写后的类型名。
std::string normalizeModelType(std::string type);

/// @brief 规整相名。
/// @param name 原始相名。
/// @return 去掉空格/下划线/连字符并转小写后的相名。
std::string normalizePhaseName(std::string name);

/// @brief 解析相角色。
/// @param role 用户配置的角色名。
/// @return 对应 PhaseRole, 无法识别时返回 Unknown。
PhaseRole parsePhaseRole(std::string role);

/// @brief 返回相角色的文本名。
std::string phaseRoleName(PhaseRole role);

/// @brief 判断模型类型是否为 Level Set。
bool isLevelSetType(const std::string& type);

/// @brief 判断模型类型是否为 mixture-like 单速度 alpha 模型。
bool isMixtureType(const std::string& type);
/// @brief 判断是否为单速度、单压、单温的均质多相 FluidStateModel。
bool isHomogeneousType(const std::string& type);
/// @brief 判断是否为每相独立速度、共享压力的 Eulerian-Eulerian 系统。
bool isEulerianEulerianType(const std::string& type);

/// @brief 判断模型类型是否为单温度热输运模型。
bool isThermalType(const std::string& type);

/// @brief 解析多相求解器家族。
MultiPhaseSolverKind parseMultiPhaseSolverKind(std::string solver);

/// @brief 返回多相求解器家族名。
std::string multiPhaseSolverKindName(MultiPhaseSolverKind solver);

/// @brief 根据 auto/模型内容解析实际求解器家族。
MultiPhaseSolverKind resolvedMultiPhaseSolverKind(
    const MultiPhaseConfig& config);

/// @brief 构造 density-based mixture 的扩展守恒量布局说明。
std::vector<std::string> densityBasedConservedVariables(
    const MultiPhaseConfig& config);

/// @brief 查找已声明相。
/// @param config 多相配置。
/// @param name 相名。
/// @return 找到时返回指针, 否则返回 nullptr。
const PhaseProperties* findPhase(const MultiPhaseConfig& config,
                                 const std::string& name);

/// @brief 在给定温度计算物性定律；常数定律返回 fallback。
double temperaturePropertyValue(const TemperaturePropertyLaw& law,
                                double fallback,
                                double temperature,
                                const std::string& propertyName);

/// @brief 计算相密度；rho(T) 优先，否则使用显式 EOS。
double phaseDensity(const PhaseProperties& phase,
                    double temperature,
                    double pressure);

/// @brief 计算动力黏度。
double phaseViscosity(const PhaseProperties& phase, double temperature);

/// @brief 计算定压比热。
double phaseSpecificHeat(const PhaseProperties& phase, double temperature);

/// @brief 计算导热系数。
double phaseThermalConductivity(const PhaseProperties& phase,
                                double temperature);

/// @brief 由 Cp(T) 积分得到比焓；常 Cp 保持旧式 h=Cp*T 基准。
double phaseSpecificEnthalpy(const PhaseProperties& phase,
                             double temperature);

/// @brief 由比焓反求温度；多项式定律仅在声明温区内求根。
double phaseTemperatureFromEnthalpy(const PhaseProperties& phase,
                                    double enthalpy);

/// @brief 返回相对的代表直径，pair 显式值优先于分散相设置。
double phasePairDiameter(const MultiPhaseConfig& config,
                         const PhasePairModelOptions& pair);

/// @brief 查找无序相对的表面张力。
double phasePairSurfaceTension(const MultiPhaseConfig& config,
                               const std::string& first,
                               const std::string& second);

/// @brief 根据相名解析 Level Set 符号。
/// @param config 多相配置。
/// @param phaseName 相名。
/// @return liquid/water 为 +1, gas/air/vapor 为 -1。
double phaseSign(const MultiPhaseConfig& config,
                 const std::string& phaseName);

/// @brief 校验多相配置, 发现不明确或不支持设置时抛出异常。
/// @param config 多相配置。
/// @param context 错误消息中的配置来源说明。
void validateMultiPhaseConfig(const MultiPhaseConfig& config,
                              const std::string& context);

} // namespace Multiphase
} // namespace Physics
} // namespace SF
