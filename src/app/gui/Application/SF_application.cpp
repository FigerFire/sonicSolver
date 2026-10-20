/// @file SF_application.cpp
/// @brief GUI 应用生命周期与主窗口启动职责。

#include "SF_application.h"

#include "SF_controller.h"
#include "SF_mainWindow.h"
#include "SF_splash.h"
// 其他 GUI 相关头文件
#include <QApplication>
// 文件系统与字体相关头文件
#include <QDir>
#include <QFileInfo>
#include <QFont>
// 设置 VTK/OpenGL 格式，以确保 VTK 渲染窗口能够正确使用 OpenGL 上下文
#include <QSurfaceFormat>
// 延迟执行任务，以便在 GUI 初始化完成后执行某些操作
#include <QTimer>

    // 判断是否启动VTK
#ifdef SONIC_GUI_WITH_VTK
#include <QVTKOpenGLNativeWidget.h>
#endif

int SF::GUI::Application::run(int argc, char* argv[])
{
// 设置默认的 OpenGL 格式，以确保 VTK 渲染窗口能够正确使用 OpenGL 上下文
#ifdef SONIC_GUI_WITH_VTK
    QSurfaceFormat vtkFormat = QVTKOpenGLNativeWidget::defaultFormat();
    // macOS 上多重采样与共面 actor 叠加时，交互旋转可能留下未清理的碎片。
    // VTK 自己负责抗锯齿；这里固定双缓冲和完整深度/模板缓冲。
    vtkFormat.setSamples(0);
    vtkFormat.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
    vtkFormat.setDepthBufferSize(24);
    vtkFormat.setStencilBufferSize(8);
    QSurfaceFormat::setDefaultFormat(vtkFormat);
#endif

    // 初始化核心应用程序
    QApplication application(argc, argv);
    // 设置应用程序的名称、组织名称和默认字体
    application.setApplicationName(QStringLiteral("sonicGui"));
    application.setOrganizationName(QStringLiteral("sonicSolver"));
    #ifdef Q_OS_MAC
    application.setFont(QFont(QStringLiteral(".AppleSystemUIFont"), 12));
    #else
    application.setFont(QFont(QStringLiteral("Segoe UI"), 10));
    #endif

    // 创建并显示启动窗口
    SplashScreen splash;
    splash.showCentered();// 显示启动窗口并居中
    splash.setStage(QStringLiteral("正在加载 Qt 与 VTK 显示模块..."));// 设置启动阶段文字1

    // 创建主窗口和 GUI 控制器
    MainWindow window;
    splash.setStage(QStringLiteral("正在装配 GUI 与求解器控制器..."));// 设置启动阶段文字2
    GuiController controller(&window, QStringLiteral(SONIC_GUI_SOURCE_DIR));// 创建 GUI 控制器，传入主窗口和源代码目录
    // SONIC_GUI_SOURCE_DIR 是 CMake 里定义的宏
    controller.initialize();// 初始化控制器
    splash.setStage(QStringLiteral("工作台已就绪"));

    // 延迟关闭启动画面并显示主窗口
    QTimer::singleShot(850, &application, [&splash, &window] {
        splash.finish(&window);// 淡出启动窗口，并激活主窗口
    });

    // 测试专用：如果设置了环境变量 SONIC_GUI_SCREENSHOT
    // 则在启动后 1.2 秒截取主窗口的屏幕截图并保存到指定路径，然后退出应用程序
    const QString screenshotPath =
        qEnvironmentVariable("SONIC_GUI_SCREENSHOT");
    if (!screenshotPath.isEmpty()) {
        QTimer::singleShot(1200, &application, [&window, screenshotPath] {
            QDir().mkpath(QFileInfo(screenshotPath).absolutePath());
            window.grab().save(screenshotPath);
            QApplication::quit();
        });
    }
    // 进入应用程序的主事件循环
    return application.exec();
}
