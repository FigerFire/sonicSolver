/// @file SF_widget.h
/// @brief GUI 内嵌 VTK 结果视图实现。

#pragma once

#include "SF_viewerTypes.h"

#include <QWidget>
#include <QVector>

#include <functional>

class QEvent;
class QDoubleSpinBox;
class QFileSystemWatcher;
class QLabel;
class QTimer;

namespace SF::GUI {

class ViewerToolbar;

/// @brief SFM/VTS/VTM/VTU/STL 网格与 IBM 几何预览部件。
///
/// 这个类是 GUI 右上角 3D 预览窗口的总入口：
/// - 读取输入网格 `.sfm`、求解器输出的 `.pvd/.vts/.vtu/.vtm` 或 STL 文件；
/// - 把不同 VTK 数据统一转换为 `vtkPolyData`，交给同一套 mapper/actor 渲染；
/// - 管理工具栏里的物理量选择、Surface/Surface With Edges、透视开关和视角按钮；
/// - 在结果目录变化时自动刷新最新时间步。
class VTKView : public QWidget {
    Q_OBJECT

public:
    /// @brief 创建 3D 预览窗口、工具栏、切片输入面板和 VTK 渲染器。
    explicit VTKView(QWidget* parent = nullptr);
    /// @brief 使用主窗口外部提供的统一工具栏创建 3D 预览窗口。
    explicit VTKView(ViewerToolbar* toolbar, QWidget* parent = nullptr);

    /// @brief 清理全局事件过滤器和 VTK 实现对象。
    ~VTKView() override;

    /// @brief 读取单个 SFM/VTK/STL 文件，或转交给 PVD 时间序列读取。
    /// @param path 文件路径，支持 `.sfm/.pvd/.vts/.vtu/.vtm/.stl`。
    /// @param error 失败时写入用户可读错误信息。
    /// @return 成功读取并刷新预览时返回 true。
    bool loadFile(const QString& path, QString* error = nullptr);

    /// @brief 将 OpenFOAM blockMeshDict 直接渲染为参数化三维几何预览。
    /// @details 这是 GeometryModel 的派生 VTK 表示，不读取或生成真实 polyMesh。
    bool loadBlockMeshGeometry(const QString& path, QString* error = nullptr);

    /// @brief 当前构建是否启用了 VTK 预览。
    bool vtkAvailable() const;

    /// @brief 监听结果目录中的 PVD，并自动显示最新时间步。
    /// @param directory case 的结果目录。
    void watchResultDirectory(const QString& directory);

    /// @brief 立即重新读取当前结果目录中的 PVD。
    void refreshResult();
    /// @brief 进入 Mesh 视图时同步默认相机和投影按钮状态。
    void activateMeshDefaults();
    /// @brief 切换为参数化 Geometry 场景，不与结果 actor 共面显示。
    void showGeometryScene();
    /// @brief 切换为 Mesh/Result 场景，不与 Geometry actor 共面显示。
    void showResultScene();
    /// @brief 注册实际显示文件变化回调，由 Controller 同步输出任务高亮。
    void setVisualizationLoadedHandler(
        std::function<void(const QString&)> handler);

signals:
    void statusMessage(const QString& message);

protected:
    /// @brief 拦截空鼠标移动，避免鼠标划过 3D 窗口时触发相机旋转。
    bool eventFilter(QObject* watched, QEvent* event) override;

public:
    void highlightBlock(int blockId);
    void highlightPatch(int patchId);
    void highlightFace(int patchId, int faceId);
    void clearHighlight();

private:
    struct TimeFrame {
        QString filePath;
        QString physicalTime;
    };

    /// @brief 读取单个真实数据文件，可选择是否重置相机。
    bool loadDataFile(const QString& path,
                      QString* error,
                      bool resetView,
                      bool fromTimeSeries);
    /// @brief 读取 PVD 时间序列并显示合适的当前时间步。
    bool loadPvd(const QString& path, QString* error);
    /// @brief 读取指定 PVD 时间步。
    bool loadTimeFrame(int frameIndex, QString* error, bool resetView);

    /// @brief 延迟刷新结果，合并 QFileSystemWatcher 的连续触发。
    void scheduleResultRefresh();

    /// @brief 恢复固定默认视角。
    void resetCamera();

    /// @brief 打开或关闭切片输入面板。
    void setSliceEnabled(bool enabled);

    /// @brief 从切片面板读取坐标，并按当前相机方向生成切片。
    void applySliceFromPanel();

    /// @brief 关闭切片并恢复完整 surface 数据显示。
    void closeSlicePanel();

    /// @brief 切换用于着色的点标量数组。
    void setScalarField(const QString& scalarName);

    /// @brief 切换物理面或 Surface With Edges 显示模式。
    void setDisplayMode(ViewerDisplayMode mode);

    /// @brief 切换透视/正交投影。
    void setPerspectiveEnabled(bool enabled);

    /// @brief 把相机切到指定轴向的正视角。
    void setCameraDirection(ViewerCameraDirection direction);
    /// @brief 开始或停止 PVD 时间步连续播放。
    void setPlaybackEnabled(bool enabled);
    /// @brief 设置连续播放速度，范围 1–999 frame/s。
    void setPlaybackRate(int framesPerSecond);
    /// @brief 播放定时器推进一帧。
    void advancePlayback();
    /// @brief 切到上一时间步。
    void showPreviousTimeStep();
    /// @brief 切到下一时间步。
    void showNextTimeStep();
    /// @brief 切到起始时间。
    void showFirstTimeStep();
    /// @brief 切到最后时间。
    void showLastTimeStep();
    /// @brief 同步工具栏时间步状态。
    void updateTimeControls();

    /// @brief 根据当前切片、物理量和显示模式重建 mapper/actor 状态。
    void updateVisualizationPipeline();

    /// @brief 从当前 surface 数据收集可选标量名，并同步到工具栏下拉框。
    void updateScalarChoices();

    /// @brief 按当前数据包围盒更新切片坐标输入范围和步长。
    void updateSlicePanelControls();

    /// @brief 使用当前切片原点和法向重新计算 vtkCutter 输出。
    void updateSliceData();

#ifdef SONIC_GUI_WITH_VTK
    /// @brief PIMPL：隔离 VTK 头文件和智能指针，降低头文件传播成本。
    class Impl;
    Impl* impl_ = nullptr;
#else
    /// @brief 未启用 VTK 时显示的占位提示。
    QLabel* placeholder_ = nullptr;
#endif
    /// @brief 3D 预览工具栏，负责发出视图控制信号。
    ViewerToolbar* toolbar_ = nullptr;
    /// @brief VTK canvas 与切片面板叠放的父容器。
    QWidget* viewerArea_ = nullptr;
    /// @brief 点击 Slice 后显示在 3D 窗口左上角的坐标输入面板。
    QWidget* slicePanel_ = nullptr;
    QDoubleSpinBox* sliceXInput_ = nullptr;
    QDoubleSpinBox* sliceYInput_ = nullptr;
    QDoubleSpinBox* sliceZInput_ = nullptr;
    /// @brief 监听结果目录和 PVD 文件变化。
    QFileSystemWatcher* resultWatcher_ = nullptr;
    /// @brief 合并文件系统事件，避免同一输出写入过程触发多次读取。
    QTimer* resultReloadTimer_ = nullptr;
    /// @brief 时间序列播放定时器。
    QTimer* playbackTimer_ = nullptr;
    QString resultDirectory_;
    QString preferredPvd_;
    QString watchedPvd_;
    QString loadedPvdPath_;
    /// @brief 当前可复用静态拓扑属于哪个 PVD；切换时间序列时强制重建一次。
    QString topologyPvdPath_;
    QString loadedDataPath_;
    QVector<TimeFrame> timeFrames_;
    int currentFrameIndex_ = -1;
    /// @brief 默认优先显示压力；没有 Pressure 时自动选择第一个可用标量。
    QString selectedScalar_ = QStringLiteral("Pressure");
    /// @brief 当前切片面经过的点。
    double sliceOrigin_[3] = {0.0, 0.0, 0.0};
    /// @brief 当前切片面的法向，应用切片时取当前相机朝向。
    double sliceNormal_[3] = {1.0, 0.0, 0.0};
    /// @brief 当前数据集包围盒，用于限制切片输入范围。
    double sliceBounds_[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    ViewerDisplayMode displayMode_ = ViewerDisplayMode::PhysicalSurface;
    bool sliceEnabled_ = false;
    /// @brief 默认使用透视投影，进入 Mesh 页时与工具栏按钮保持一致。
    bool perspectiveEnabled_ = true;
    int playbackRate_ = 24;
    /// @brief 用于判断 PVD 引用的数据文件是否已经变化。
    qint64 loadedDataModified_ = -1;
    qint64 loadedDataSize_ = -1;
    std::function<void(const QString&)> visualizationLoadedHandler_;
};

} // namespace SF::GUI
