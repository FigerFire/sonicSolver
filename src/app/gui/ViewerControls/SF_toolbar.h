#pragma once

/// @file SF_toolbar.h
/// @brief 3D 可视化窗口工具栏。

#include "SF_viewerTypes.h"

#include <QList>
#include <QWidget>

#include <functional>

class QComboBox;
class QHBoxLayout;
class QLabel;
class QSlider;
class QWidget;
class QToolButton;

namespace SF::GUI {

/// @brief 类 ParaView 的 3D 预览工具栏。
///
/// 工具栏只负责 UI 状态和信号，不直接接触 VTK 对象。
class ViewerToolbar final : public QWidget {
    Q_OBJECT

public:
    explicit ViewerToolbar(QWidget* parent = nullptr);

    /// @brief 根据当前构建是否支持 VTK 启停工具按钮。
    void setVtkAvailable(bool available);
    /// @brief 把主窗口的模块下拉和 GEO/MESH 开关嵌入第一行工具栏。
    void setNavigationControls(QWidget* moduleSelector,
                               QWidget* viewModeSwitch);
    /// @brief 在 GEO/MESH 两套工具组之间切换。
    void setWorkspaceMode(ViewerWorkspaceMode mode);
    /// @brief 更新可显示的点数据标量列表。
    /// @param scalars 标量数组名。
    /// @param preferredScalar 优先选择的数组名，通常是 Pressure。
    void setScalarChoices(const QStringList& scalars,
                          const QString& preferredScalar);
    /// @brief 同步 slice 按钮状态。
    void setSliceEnabled(bool enabled);
    /// @brief 同步显示模式按钮状态。
    void setDisplayMode(ViewerDisplayMode mode);
    /// @brief 同步投影按钮状态；true 表示透视投影。
    void setPerspectiveEnabled(bool enabled);
    /// @brief 同步六向视角按钮状态。
    void setViewDirection(ViewerCameraDirection direction);
    /// @brief 同步时间步显示和可用按钮。
    void setTimeFrameState(int currentIndex,
                           int frameCount,
                           const QString& physicalTime);
    /// @brief 同步播放按钮状态。
    void setPlaybackEnabled(bool enabled);
    /// @brief 注册播放帧率变化回调。
    void setPlaybackRateHandler(std::function<void(int)> handler);

signals:
    /// @brief 用户切换 slice 显示。
    void sliceToggled(bool enabled);
    /// @brief 用户选择新的物理标量。
    void scalarChanged(const QString& scalarName);
    /// @brief 用户切换“物理面/网格”显示模式。
    void displayModeChanged(ViewerDisplayMode mode);
    /// @brief 用户切换透视/正交投影。
    void perspectiveChanged(bool enabled);
    /// @brief 用户点击六向视角按钮。
    void viewDirectionRequested(ViewerCameraDirection direction);
    /// @brief 用户请求重置相机到适应窗口。
    void resetCameraRequested();
    /// @brief 用户切换时间序列连续播放。
    void playbackToggled(bool enabled);
    /// @brief 用户请求上一时间步。
    void previousTimeStepRequested();
    /// @brief 用户请求下一时间步。
    void nextTimeStepRequested();
    /// @brief 用户请求回到起始时间。
    void firstTimeStepRequested();
    /// @brief 用户请求跳到最后时间。
    void lastTimeStepRequested();
    /// @brief 用户选择几何草图工具；当前只负责 UI 和状态。
    void geometryToolChanged(ViewerGeometryTool tool);
    /// @brief 用户输入拉伸长度。
    void solidExtrudeRequested(double length);
    /// @brief 用户输入旋转扫描角度。
    void solidRevolveRequested(double angleDegrees);

private:
    QToolButton* makeIconButton(const QString& iconPath,
                                const QString& tooltip,
                                bool checkable = false);
    QToolButton* makeToolButton(const QString& iconPath,
                                const QString& fallbackText,
                                const QString& tooltip,
                                bool checkable = false);
    void addViewButton(ViewerCameraDirection direction,
                       const QString& iconPath,
                       const QString& tooltip);
    void addGeometryTool(ViewerGeometryTool tool,
                         const QString& iconPath,
                         const QString& tooltip);
    void registerVtkControl(QWidget* control);
    void registerMeshControl(QWidget* control);
    void registerGeometryControl(QWidget* control);

    QHBoxLayout* firstRowLayout_ = nullptr;
    QHBoxLayout* secondRowLayout_ = nullptr;
    QHBoxLayout* thirdRowLayout_ = nullptr;
    QList<QWidget*> vtkControls_;
    QList<QWidget*> meshControls_;
    QList<QWidget*> geometryControls_;
    QToolButton* sliceButton_ = nullptr;
    QToolButton* surfaceButton_ = nullptr;
    QToolButton* meshButton_ = nullptr;
    QToolButton* perspectiveOnButton_ = nullptr;
    QToolButton* perspectiveOffButton_ = nullptr;
    QToolButton* resetButton_ = nullptr;
    QToolButton* playButton_ = nullptr;
    QToolButton* previousTimeButton_ = nullptr;
    QToolButton* nextTimeButton_ = nullptr;
    QToolButton* firstTimeButton_ = nullptr;
    QToolButton* lastTimeButton_ = nullptr;
    QToolButton* extrudeButton_ = nullptr;
    QToolButton* revolveButton_ = nullptr;
    QLabel* timeLabel_ = nullptr;
    QLabel* playbackRateLabel_ = nullptr;
    QSlider* playbackRateSlider_ = nullptr;
    QComboBox* scalarSelector_ = nullptr;
    std::function<void(int)> playbackRateHandler_;
};

} // namespace SF::GUI
