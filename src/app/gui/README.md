# sonicGui architecture

`sonicGui` 使用分层静态库组织，界面层不直接解析求解器内部状态。

- `GUI/`: 窗口、按钮、输入框、配置树和参数编辑控件。
- `hotKey/`: 单键热键路由，把 G/M/S 映射到几何、网格、求解器模块。
- `Geometry/`: 单块二维截面的线段/圆弧草图和网格参数生成。
- `Controller/`: 接收界面信号，协调配置 IO、求解器进程、终端和预览。
- `ViewerControls/`: 3D 预览工具栏、显示模式枚举和 ParaView 风格鼠标交互。
- `VTKView/`: VTS、VTU、VTM、STL 的显示、线框切换和相机控制。
- `IO/`: 工作项目和 OpenFOAM-like 文档模型，负责 `system/` 字典、模块模板与导入文件。
- `Application/`: `SF_application.*` 管理 Qt 应用生命周期，`main.cpp` 只转发。
- `Resources/`: 普通 UI 图标、配色素材和资源入口说明；按钮、菜单、树节点图标都从这里走。
- `Splash/assets/frames/`: 启动动画帧，只放 VTK 压力切片这类启动图，不放普通按钮图标。

CFD 求解和正式 case/mesh/result IO 不在 GUI 内重复实现：

- `src/solver/algorithm/`: 显式/隐式、密度基/压力基的求解循环。
- `src/solver/equation/`: 方程定义和装配。
- `src/solver/discretization/`、`src/solver/boundary/`: 离散和边界策略。
- `src/app/application/`: case 装配和运行入口。
- `src/infrastructure/io/`: 求解器使用的 case、mesh、SFM、VTK/VTM IO。

`src/app/application/SF_application.*` 不包含 `main()`。
命令行入口位于 `src/app/cli/main.cpp`，GUI 入口位于
`src/app/gui/Application/main.cpp`。GUI 通过 `QProcess` 启动命令行前端，
以隔离 MPI 生命周期、计算崩溃和停止操作。

用户只选择工作文件夹。GUI 直接以 OpenFOAM-like case 目录为项目根；
左侧固定显示几何、网格、求解器、湍流模型、初始条件、边界条件和输出，
IBM、MRF、Parallel 通过模块菜单或树控件右键新增。

## 模块代码索引

`src/app/gui` 下每个目录都按“调度入口 + 实现文件”的方式组织：

- `Application/`
  - `SF_application.h/.cpp`: GUI 应用生命周期，创建 `QApplication`、Splash、
    `MainWindow` 和 `GuiController`。
  - `main.cpp`: GUI 可执行程序入口，只转发到 `Application::run()`。
- `Controller/`
  - `SF_controller.h/.cpp`: 界面信号到项目 IO、求解器进程、终端命令、
    网格生成和结果刷新的协调层。
- `GUI/`
  - `SF_GUI.h`: 控件层调度入口。
  - `SF_GUI.cpp`: 主窗口构造、菜单栏、信号连接、模块图标、全局样式。
  - `SF_mainWindow.h`: `MainWindow` 对外接口、信号、内部控件成员。
  - `SF_mainWindowLayout.cpp`: 左右布局、几何/网格视图栈、底部输出/终端/日志。
  - `SF_mainWindowTree.cpp`: 左侧配置树、模块展开、双击编辑、右键新增、已启用菜单。
  - `SF_switchButton.*`: true/false 配置项用的开关控件。
  - `SF_viewModeSwitch.*`: 顶部 GEO/MESH 双态切换控件。
- `hotKey/`
  - `SF_hotKey.h/.cpp`: 主窗口热键路由。当前注册 `G` 几何、`M` 网格、
    `S` 求解器；IBM/MRF 暂不注册。热键命中后复用
    `MainWindow::focusModule(...)`，效果等同顶部模块下拉选择。
- `ViewerControls/`
  - `SF_viewerControls.h`: 3D 预览控制模块调度入口。
  - `SF_toolbar.h/.cpp`: 类 ParaView 工具栏，包含 Slice、物理数据选择、
    物理面/网格显示、透视/正交投影和六向视图按钮。
  - `SF_interactorStyle.h/.cpp`: VTK 鼠标交互样式，左键拖拽旋转、右键拖拽平移。
  - `SF_viewerTypes.h`: 工具栏和 VTKView 共享的显示模式、视角方向枚举。
- `Geometry/`
  - `SF_Geometry.h`: 几何模块调度入口。
  - `SF_editor.h/.cpp`: 几何草图编辑器，负责点、边、圆弧和草图保存请求。
- `IO/`
  - `SF_IO.h`: GUI 配置 IO 调度入口。
  - `SF_configDocument.h/.cpp`: GUI 配置文档解析、条目修改和写回。
  - `SF_project.h/.cpp`: OpenFOAM-like 工作目录、模块文件、导入文件和模板管理。
- `VTKView/`
  - `SF_VTKView.h`: VTKView 调度入口。
  - `SF_widget.h/.cpp`: VTS/VTU/VTM/STL 预览、线框开关、结果目录监听和刷新。
- `Resources/`
  - `SF_resources.qrc`: GUI 图标资源注册。
  - `icons/`: 模块菜单、树节点、按钮和窗口图标。
- `Splash/`
  - `SF_splash.h/.cpp/.qrc`: 启动画面与启动动画资源。

## 资源与图标

当前 GUI 的模块下拉、左侧树节点和工作目录按钮图标都从
`Resources/SF_resources.qrc` 打包。你自己画好的图标统一放在
`Resources/icons/`，然后在 qrc 中通过 alias 暴露为
`:/sonic/icons/<name>.png` 或 `:/sonic/icons/<name>.svg`。

建议优先使用 SVG。SVG 不依赖 DPI，放大缩小最干净，适合左侧树、工具栏、
下拉菜单和按钮图标。若必须用 PNG，按显示尺寸准备多倍图：

- 小树节点/菜单图标：16px 显示，准备 `16x16`、`32x32(@2x)`。
- 工具栏/按钮图标：20-24px 显示，准备 `24x24`、`48x48(@2x)`。
- 启动画面或大卡片图片：按实际控件像素导出，当前 splash 动画帧单独放在
  `Splash/assets/frames/`，按 `frame-00.png` 这种顺序名替换即可。

PNG 文件的“DPI 数字”本身不重要，Qt 更关心像素尺寸和设备像素比。macOS
Retina 屏优先准备 2x 资源；例如界面上显示 24px，就导出 48x48 的高清图。
