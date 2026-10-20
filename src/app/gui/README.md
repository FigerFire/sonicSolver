# sonicGui native architecture

GUI 只操作已注册的 SonicFile native case；不再创建或读取 OpenFOAM-like 配置、dictionary、Allrun/backend 或 .foam。数学与执行 authority 来自同一个 WHAT/STATE/HOW/WHICH compiler。

|目录|真实职责|
|---|---|
|Application|Qt 应用、Splash、MainWindow/Controller 生命周期|
|GUI|窗口、注册文档树、原生文本编辑、configured-enabled 菜单、结果任务与进度|
|IO|CaseProject/ConfigDocument：CaseIO reader/writer、共享 native creation、staged validation/rollback|
|Controller|typed case validation、QProcess 的 native serial/MPI、一步预览、mesh生成和结果加载|
|Geometry|点/线/arc草图与native typed参数视图；不读取blockMeshDict|
|VTKView|SFM/STL/VTS/VTU/VTM/PVD；Pressure/native p默认着色、相机/切片/时间帧|
|ViewerControls/hotKey|显示控制与输入路由，不解释CFD数学|
|Resources/Splash|图标、配色、启动动画；无runtime authority|

`CaseProject::open` 对空目录调用与 CLI init 相同的 `Application::createNativeCase`。非空旧格式目录拒绝；不覆盖已有项目。默认只有可信 compressibleSingleFluid 模板，其余 recipe Unsupported。文件 object 与目录对应，solvers 目录的 object 是 solver；所有子文件按唯一 registry显式注册，child identity与role双重验证。

导航来自 CaseIO 真正加载的文档、文件头identity和typed对象。未注册文件不加载；模型存在与enabled:false区分。新增optional turbulence/IBM/MRF默认inactive并更新registry。复杂WHAT/HOW采用原生YAML文本编辑，保存先经过生产 reader/ModelLoader/inspection，再原子替换，失败恢复原文件；不使用regex模拟YAML AST，不新增GUI数学注册表。删除模块可编辑registry，保存仍检查其引用与数学依赖。

网格参数从registered generation（structuredMesh role）取得；同一MeshParameters冻结后由MeshGen生成SFM。GUI四角挤出/arc草图有native参数实现，复杂拓扑编辑UI明确Unsupported。SFM导入使用生产mesh loader验证，普通import不覆盖已有payload；显式regeneration可以原子替换并重新注册。STL导入是geometry contribution，不自动打开IBM。原几何/TFI数学保留。

GUI运行是 `sonicSolver run CASE`，MPI参数来自同一typed runtime，经QProcess隔离。预览复制到仓库test/guiNativePreview，执行 `sonicSolver steps 1 CASE`；正常完成后加载PVD，不能以首个输出事件杀进程来伪造一步。PVD标签为提交后的physical time，step只用于文件编号；任务/进度读取真实runtime与PVD。

构建：`./sonicMake.sh -gui`。需要自定义树时支持 SONIC_BUILD_DIR、SONIC_CMAKE_GENERATOR、SONIC_JOBS；正常默认流程不变。Qt6/VTK9.6已在本轮构建验证；离屏平台不能创建OpenGL时只算Qt/QProcess验证，host Cocoa窗口截图才证明实际VTK渲染。

测试hooks仅用于仓库test/验收：SONIC_GUI_WORKSPACE、SONIC_GUI_TEST_ACTION（run/preview/mesh）、SONIC_GUI_TEST_RESULT、SONIC_GUI_SCREENSHOT。QtCore工程harness拒绝仓库test/之外的工作目录。实际source/命令/配置/结果与能力边界见[本轮报告](../../../docs/native-architecture-closure.md)及[IO契约](../../../docs/native-io-gui-contract.md)。
