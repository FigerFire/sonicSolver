# AGENTS.md — SonicSolver Development Rules

> SonicSolver 不再以“选择一个 CFD solver”为核心，而是组合一个可解释、可执行的数值系统。
>
> **核心原则：主未知量、方程、模块贡献、耦合、时间推进、空间离散、运行时并行语义必须解耦。**

---

## 1. 项目目标

SonicSolver 是科研型结构网格 FDM 框架。一个 case 应由以下独立问题组成：

```text
1. 解哪些主未知量？
2. 解哪些方程？
3. 模块给已有方程增加什么项？
4. 模块是否增加新的方程或约束？
5. 方程之间如何耦合？
6. 时间如何推进？
7. 每个算子如何离散？
8. 串行/MPI 如何保证同一个离散系统？
```

禁止把这些问题重新打包成：

```text
compressibleSolver
incompressibleSolver
densityBasedEulerianSolver
pressureBasedMultiphaseSolver
```

---

## 2. 所有扩展模块只回答三个问题

任何新模块必须先说明自己属于以下一种或多种数学作用。

### 2.1 给已有方程增加项

例如：

```text
gravity
MRF
wall heat
surface tension
Lee phase change
RPI
explicit IBM forcing
```

数学意义：

```text
已有方程 + source / flux / coefficient contribution
```

例如：

```text
ddt(rhoU)
+ div(momentumFlux)
+ diffusion(mu,U)
=
gravity
+ MRF
+ IBM
```

**源项模块不能拥有主时间循环。**

---

### 2.2 增加方程或约束

例如：

```text
turbulence transport equation
species equation
phase equation
pressure-correction equation
level-set equation
DLM/KKT constraint
```

模块必须明确声明：

```text
新增 unknown
新增 equation / constraint
依赖哪些已有 unknown
给哪些已有方程贡献 closure/source
```

禁止把“增加一个方程”直接实现成一个新的 monolithic solver family。

---

### 2.3 改变“怎么解方程”

例如：

```text
Euler / RK2 / SSPRK3 / RK4
SIMPLE / PISO / PIMPLE
Newton / monolithic / KKT
WENO3 / WENO5 / WENO7 / TENO5
Rusanov / Steger-Warming / Roe / HLLC
Central diffusion
HYPRE backend
```

这些属于数值策略，不属于物理模型 ownership。

---

## 3. 正交的架构维度

以下概念必须保持独立。

### 3.1 Primary-variable formulation

当前最重要的两类：

```text
density
pressure
```

它只回答：

> 哪一组变量是 primary numerical unknowns？

典型 density formulation：

```text
Q = [rho, rhoU, rhoE, ...]
```

典型 pressure formulation：

```text
Q = [p, U, h/T, ...]
```

禁止硬绑定：

```text
density == compressible
pressure == incompressible
```

可压缩/不可压缩由**方程和 closure**决定，不是顶层 solver 开关。

---

### 3.2 Physical / phase system

例如：

```text
SingleFluid
HomogeneousMultiphase
OneFluidInterface
EulerianEulerian
EulerianLagrangian
```

它回答：

> 系统中有哪些连续体/相/物理状态？

它不负责选择：

```text
density / pressure
RK / BDF
SIMPLE / PIMPLE
WENO / TENO
```

---

### 3.3 Coupling algorithm

例如：

```text
none
SIMPLE
PISO
PIMPLE
Newton
monolithic
KKT
```

它回答：

> 相互依赖的方程/约束如何达到一致？

`SIMPLE/PISO/PIMPLE` 不得永久绑定到 pressure formulation。

如果用户选择：

```text
formulation = density
coupling = SIMPLE
```

框架不能因为“这不常见”而拒绝。

只允许因为**实现能力不支持**而 fail fast。

---

### 3.4 Time integration

例如：

```text
Euler
RK2
SSPRK3
RK4
BDF2
```

它回答：

> 时间导数怎么推进？

它和 coupling 是不同层级。

---

### 3.5 Spatial discretization

例如：

```text
WENO3
WENO5
WENO7
TENO5
Central2
Central4
```

通量构造/分裂又是另一选择：

```text
Rusanov
Lax-Friedrichs
Steger-Warming
Roe
HLLC
```

离散模块不得从 solver 名称推断 physical system。

---

## 4. 不设置 compressible / incompressible 顶层开关

禁止新增：

```yaml
compressible: true
```

或：

```yaml
flowType: incompressible
```

物理由方程决定。

例如：

```text
ddt(rho) + div(rho U) = 0
p = p(rho,T)
```

自然描述可压缩系统。

例如：

```text
div(U) = 0
rho = constant
```

自然描述不可压缩约束。

框架不能因为这些方程自动替换用户选择的 formulation。

---

## 5. 默认用户是专家

SonicSolver 不负责阻止“物理上不聪明”的组合。

Validator 只检查：

> 当前实现能不能执行这套组合？

允许的 fail-fast：

```text
WENO characteristic reconstruction 被选择，
但 active EOS 没有提供所需 eigenstructure。

SIMPLE 被选择，
但系统没有 pressure-correction relation/equation。

RK4 被选择，
但 active algebraic constraint 没有 stage treatment。

Monolithic KKT 被选择，
但 distributed constraint row ownership 尚未实现。
```

禁止：

```text
用户选 pressure formulation -> 自动改成 PIMPLE
用户选 density formulation -> 自动启用 PerfectGas compressible path
```

---

## 6. Equation 层的职责

Equation 层只回答：

> 正在解什么数学方程？

目标形式：

```text
EquationSystem
├── Mass
├── Momentum
├── Energy
├── Phase equations
├── Turbulence equations
└── Constraints
```

Equation 不拥有全局 timestep lifecycle。

方程应该由显式 term 组成：

```text
ddt(...)
div(...)
gradient(...)
diffusion(...)
source(...)
constraint(...)
algebraic closure
```

最终同一份 equation description 应同时服务：

```text
validation
sonicSolver explain
assembly
execution
```

不要长期保留“一套只用于打印的 DSL + 一套完全独立的手写执行逻辑”。

---

## 7. Source 模块规则

Source module 只提供数学贡献。

典型语义：

```text
source.evaluate(state, targetEquation, time)
```

或项目中已有的等价机制。

例如：

```text
gravity -> momentum RHS
MRF -> momentum RHS
Lee -> phase mass + energy
RPI -> wall/phase/energy
surface tension -> momentum
```

Source module 禁止：

```text
推进全局时间
运行 RK stages
拥有 MPI barrier
拥有主 equation lifecycle
```

---

## 8. 增加方程的模块规则

增加 equation/constraint 的模块必须声明：

```text
unknowns
new equations / constraints
required fields
required closures
coupling relationships
```

例如 turbulence：

```text
adds: k, omega
adds: k-equation, omega-equation
contributes: mu_t closure
```

例如 EulerianEulerian：

```text
adds: per-phase state
adds: phase mass/momentum/energy equations
adds: interphase exchange
```

例如 DLM/KKT：

```text
adds: lambda
adds: J u = U_b
adds: J and J^T algebraic coupling
```

---

## 9. Coupling 规则

Coupling algorithm 操作的是**已经声明好的 equations**。

```text
none:
    不增加额外 consistency iteration

SIMPLE:
    predictor -> correction -> repeat/commit

PISO:
    predictor -> one/multiple correction passes

PIMPLE:
    outer nonlinear loop + PISO-like inner correction

KKT:
    primal + constraint coupled solve
```

Coupling 不得静默重定义物理方程。

---

## 10. Time Integrator 规则

Time integrator 应只知道：

```text
stage coefficients
stage times
which state is integrated
how to request RHS/equation evaluation
how to publish stage state
```

它不应该知道：

```text
WENO 内部实现
MRF 实现
IBM geometry
MPI communicator
具体 phase model 类名
```

---

## 11. Discretization 规则

Discretization 只回答：

> 一个数学 operator 怎么变成数值贡献？

例如：

```text
div(momentumFlux)
    -> WENO7
    -> characteristic reconstruction
    -> Steger-Warming
```

```text
diffusion(mu,U)
    -> Central2
```

```text
ddt(Q)
    -> RK4 stage relation / implicit matrix contribution
```

Discretization 不拥有全局时间循环。

---

## 12. IBM 语义

不同 IBM 方法必须按数学作用分类。

### Ghost / ILW

```text
boundary closure
```

执行位置：

```text
physical boundary
-> halo
-> ghost/ILW reconstruction
-> spatial operator reads ghost state
```

### Explicit / fractional forcing

根据实际算法属于：

```text
equation source contribution
```

或：

```text
post-predictor state correction
```

不能为了统一 API 假装二者完全一样。

### DLM / KKT / fully implicit constraint

```text
adds unknown(s)
adds constraint equation(s)
changes coupled algebra
```

属于：

```text
Equation System + Coupling + Linear Algebra
```

不能拥有第二套 CFD timestep lifecycle。

---

## 13. 并行语义是数值不变量

Raw MPI 只能存在于 infrastructure/backend。

### State / geometry / labels

```text
owner -> COPY -> replicas
```

禁止平均：

```text
Q
geometry
metrics
cell type
classification
ownership metadata
```

### Shared numerical face

一个物理 face 只有一个 canonical numerical flux：

```text
candidate information ready
-> one canonical F*
-> COPY F*
-> owner +F*
-> neighbour -F*
```

禁止两边独立算 F* 再平均。

### Residual / source / load / constraint contribution

```text
local contributors
-> SUM
-> canonical owner
```

禁止把 SUM 改成 COPY/average。

---

## 14. State 与 Workspace ownership

每一种语义对象只有一个 authority。

### Physical state

禁止引入：

```text
第二份 Q authority
重复 clocks
重复 equation bindings
重复 phase registries
```

### Solver workspace

例如：

```text
FluxField
Residual
RK k1/k2/k3/k4
linear-system workspace
constraint solve workspace
```

属于 solver execution lifetime。

不要重新塞回 `Field`。

---

## 15. 当前源码的过渡映射

当前代码仍处于迁移期，可以暂时这样理解：

```text
CompressibleAlgorithm
    -> 当前 density-formulation driver
    -> 长期不应与“compressible”绑定

DensityBasedTime
    -> 当前 explicit time integration

DensityBasedRHS
    -> 当前 explicit equation assembly executor

Equation::Compressible::System
    -> 当前 mass/momentum/energy equation bundle

EulerianEulerian::PressureStepper
    -> 当前混合了：
       physical system
       primary formulation
       coupling
       equation execution
    -> 是长期拆分对象

PatchWorkspace
    -> solver-owned numerical workspace

StateBundle
    -> state composition / clock / equation binding

ResolvedSimulationSystem
    -> 当前 description/validation representation
    -> 优先演化为 authoritative executable/explainable system
```

除非证明现有结构无法承载，否则不要新建平行的 `SimulationPlan` 大体系。

---

## 16. 依赖主干

目标依赖方向：

```text
Application / Composition
        ↓
Resolved Numerical System
        ↓
Algorithm / Coupling / Time
        ↓
Equation
        ↓
Discretization + Boundary
        ↓
Linear Algebra
        ↓
Runtime interfaces
        ↓
Infrastructure backend
```

Models 只提供：

```text
closures
coefficients
source terms
new equations / constraints
```

Models 不拥有全局 solver lifecycle。

---

## 17. `methods/` 规则

`methods/` 只放 solver-independent math / geometry kernels。

合适：

```text
polynomial math
small matrix math
geometry utility
generic interpolation
generic root finding
```

CFD-specific execution 最终应属于 solver：

```text
WENO/TENO convection
scalar PDE transport lifecycle
viscous PDE assembly
time integration
```

但不要为了目录漂亮而制造 reverse dependency。

---

## 18. IO 原则

配置应表达独立概念，而不是 solver identity。

目标：

```yaml
formulation:
  type: density

system:
  type: eulerianEulerian

coupling:
  type: none

time:
  type: RK4

numerics:
  convection:
    reconstruction: WENO7
    flux: StegerWarming
```

物理模型属于 models：

```yaml
transport:
  type: constant

gravity:
  enabled: true

MRF:
  enabled: false
```

---

## 19. Optional Model enable 统一规则

```text
block absent
    -> disabled

block present
    -> enabled by default

block present + enabled: false
    -> configured but inactive
```

统一使用：

```yaml
enabled: true|false
```

不要混用：

```text
active
switch
on
enable
```

不要在 cell loop 里每次判断 enabled；inactive model 应在 composition 阶段退出执行路径。

---

## 20. Built-in Field Registration

模型/方程自己声明 required field。

例如：

```text
OneFluidInterface -> phi
kOmegaSST -> k, omega
EulerianEulerian -> phase states
```

用户只提供：

```text
initialization
boundary condition
model parameters
```

用户自定义 field 才显式注册。

禁止按字段名字写硬编码 auto-register：

```cpp
if (name == "phi") ...
```

---

## 21. `sonicSolver explain` 是架构约束

任何 executable configuration 都必须可以由**同一份 authoritative system description**逐层解释。

### Level 1 — Global solution flow

例如：

```text
Q^n
 ↓
assemble R(Q)
 ↓
dQ/dt = -R
 ↓
RK4
 ↓
Q^(n+1)
```

### Level 2 — Active equations

```text
Mass:
ddt(rho) + div(massFlux) = 0

Momentum:
ddt(rhoU)
+ div(momentumFlux)
+ diffusion(mu,U)
= gravity + MRF + IBM

Energy:
ddt(rhoE)
+ div(energyFlux)
+ diffusion(k,T)
= energySource
```

### Level 3 — Added equations / constraints / coupling

```text
Coupling: none / SIMPLE / PISO / PIMPLE
Constraint: DLM/KKT
```

### Level 4 — Time integration

```text
Euler / RK2 / SSPRK3 / RK4 / ...
```

### Level 5 — Operator discretization

```text
div(momentumFlux)
  reconstruction : WENO7
  variable space : characteristic
  flux split     : StegerWarming
```

### Level 6 — Boundary / parallel assembly

```text
Ghost IBM:
physical BC -> halo -> ILW -> operator

Shared face:
candidate -> canonical F* -> COPY -> +/- assembly

GlobalDof:
local residual -> SUM -> owner
```

如果 Explain 不能按这个顺序说明 case，说明 runtime architecture 仍然在混层。

---

## 22. 禁止 hidden fallback

禁止静默替换：

```text
WENO7 -> Rusanov
StegerWarming -> LaxFriedrichs
PIMPLE -> SIMPLE
density -> pressure
generic EOS -> PerfectGas
```

不支持就 fail fast，并明确说明缺少什么 capability。

---

## 23. 数值重构验收

结构重构不能因为：

```text
build 成功
不 crash
代码更漂亮
```

就算完成。

必须保护 frozen numerical baseline。

数值路径修改至少比较：

```text
dt
stage times
mass
momentum
energy
rho/p extrema
norms
canonical face behavior
GlobalDof behavior
```

出现 regression 时寻找**第一个不同的 checkpoint**：

```text
initial Q
dt
boundary-ready Q
candidate flux
canonical flux
directional residual
local residual
global residual
stage Q
```

---

## 24. 重构优先级

默认：

```text
delete
>
simplify
>
reuse
>
add abstraction
```

不要为了“架构感”新建无必要的：

```text
Manager
Adapter
Context
Registry
Service
Descriptor
Coordinator
Facade
Strategy
```

只有当它代表真实数学概念、消除真实依赖、或支持多个真实实现时才允许新增抽象。

---

## 25. 注释规则

注释只解释：

```text
WHY
MATHEMATICAL MEANING
OWNERSHIP
LIFECYCLE
PARALLEL INVARIANT
```

好：

```cpp
// One physical interface owns one canonical F*.
// Replicas copy it; they never average independent flux evaluations.
```

差：

```cpp
// Loop over faces.
```

---

## 26. 每次非平凡修改必须汇报

```text
1. 影响哪一层
2. 属于“加项 / 加方程 / 改求解方式”中的哪一类
3. ownership before/after
4. execution order before/after
5. public API before/after
6. dependency before/after
7. MPI/parallel semantics 是否变化
8. numerical semantics 是否变化
9. regression tests
10. deferred issues
```

---

## 27. Agent 编码前必须回答

```text
这是哪一层？

它是：
1. 加一个 term？
2. 加 equation/constraint？
3. 改怎么解 equation？

当前调用链是什么？
谁拥有 state？
谁拥有 workspace？
是否改变 stage timing / flux / residual / MPI sync？
现有架构能否直接表达，还是确实需要新 abstraction？
```

---

## 28. Hard Prohibitions

禁止：

```text
重新出现 single/multi 两套 density lifecycle
把 EulerianEulerian 永久绑定 pressure formulation
把 SIMPLE/PISO/PIMPLE 永久绑定 pressure formulation
把 RK 永久绑定 density formulation
把 compressible/incompressible 变成必选 solver family
用 fallback 掩盖不支持组合
在 equation/discretization 中直接写 MPI implementation
平均 canonical state/geometry/face flux identity
把 FluxField/Residual stage workspace 塞回 Field
创建第二个 physical-state authority
继续扩展只用于打印、无法驱动 execution 的 Equation DSL metadata
```

---

## 29. 成功标准

SonicSolver 的 case 最终应被理解为：

```text
Primary variables
      +
Equations
      +
Additional term/equation/constraint contributions
      +
Coupling
      +
Time integration
      +
Discretization
      +
Runtime execution
```

而不是：

```text
choose one pre-packaged solver name
```

---

## 30. Built-in Time Recipes

时间算法采用不可变的 built-in recipe：

```text
forwardEuler
SSPRK3
classicalRK4
```

用户选择完整 recipe，不单独设置 stage 数或 explicit/implicit flag。方程层只描述
`ddt/div/diffusion/source/constraint`；`ResolvedSimulationSystem::timeRecipe`
是时间方法与 stage 数的唯一 authority，`SolvePlanner` 只按 recipe topology 和
stageCount 生成 Plan，`Time::Explicit` 只实现单 stage 数学。

```text
Equation-driven
    + Recipe-compiled
    + Plan-executed
```

禁止恢复 `TimeScheme`、`explicitStageCount`、`ddtDispatch`，也禁止创建
`ExplicitSolver`、`ImplicitSolver`、`RK4Solver` 等全局 lifecycle 类型。
