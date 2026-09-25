/// @file SF_toolbar.cpp
/// @brief 三维视图交互与工具栏控制实现。

#include "SF_toolbar.h"

#include <QButtonGroup>
#include <QComboBox>
#include <QHBoxLayout>
#include <QIcon>
#include <QInputDialog>
#include <QLabel>
#include <QSignalBlocker>
#include <QSlider>
#include <QToolButton>
#include <QVariant>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <utility>

namespace SF::GUI {
namespace {

constexpr int playbackSliderMinimum = 0;
constexpr int playbackSliderMaximum = 999;
constexpr int playbackRateMinimum = 1;
constexpr int playbackRateMaximum = 999;

int playbackRateForSliderPosition(int position)
{
    const int boundedPosition = std::clamp(
        position, playbackSliderMinimum, playbackSliderMaximum);
    const double fraction = static_cast<double>(
        boundedPosition - playbackSliderMinimum)
        / static_cast<double>(playbackSliderMaximum
                              - playbackSliderMinimum);
    const double rate = std::exp(
        std::log(static_cast<double>(playbackRateMaximum)) * fraction);
    return std::clamp(static_cast<int>(std::lround(rate)),
                      playbackRateMinimum, playbackRateMaximum);
}

int sliderPositionForPlaybackRate(int framesPerSecond)
{
    const int boundedRate = std::clamp(
        framesPerSecond, playbackRateMinimum, playbackRateMaximum);
    const double fraction = std::log(static_cast<double>(boundedRate))
        / std::log(static_cast<double>(playbackRateMaximum));
    return static_cast<int>(std::lround(
        playbackSliderMinimum
        + fraction * static_cast<double>(playbackSliderMaximum
                                         - playbackSliderMinimum)));
}

QString formatPhysicalTime(const QString& physicalTime)
{
    if (physicalTime.trimmed().isEmpty()) return QStringLiteral("-");
    bool ok = false;
    const double value = physicalTime.toDouble(&ok);
    if (!ok) return physicalTime;
    return QString::number(value, 'e', 2);
}

} // namespace

ViewerToolbar::ViewerToolbar(QWidget* parent)
    : QWidget(parent)
{
    setObjectName(QStringLiteral("viewerToolbar"));
    setAttribute(Qt::WA_StyledBackground, true);
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(10, 5, 10, 5);
    root->setSpacing(4);

    auto makeRow = [this, root]() {
        auto* row = new QWidget(this);
        auto* layout = new QHBoxLayout(row);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setSpacing(5);
        root->addWidget(row);
        return layout;
    };

    firstRowLayout_ = makeRow();
    secondRowLayout_ = makeRow();
    thirdRowLayout_ = makeRow();

    auto* displayGroup = new QButtonGroup(this);
    displayGroup->setExclusive(true);
    surfaceButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/surfaceView.png"),
        QStringLiteral("渲染物理面"), true);
    meshButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/meshView.png"),
        QStringLiteral("显示 Surface With Edges"), true);
    displayGroup->addButton(surfaceButton_);
    displayGroup->addButton(meshButton_);
    firstRowLayout_->addWidget(surfaceButton_);
    firstRowLayout_->addWidget(meshButton_);
    registerVtkControl(surfaceButton_);
    registerVtkControl(meshButton_);
    registerMeshControl(surfaceButton_);
    registerMeshControl(meshButton_);
    firstRowLayout_->addStretch();

    scalarSelector_ = new QComboBox(this);
    scalarSelector_->setObjectName(QStringLiteral("viewerScalarSelector"));
    scalarSelector_->setFixedSize(118, 28);
    scalarSelector_->setToolTip(QStringLiteral("选择显示的物理数据"));
    secondRowLayout_->addWidget(scalarSelector_);
    registerVtkControl(scalarSelector_);
    registerMeshControl(scalarSelector_);

    auto* viewGroup = new QButtonGroup(this);
    viewGroup->setExclusive(true);
    addViewButton(ViewerCameraDirection::PositiveZ,
                  QStringLiteral(":/sonic/icons/pZ.png"),
                  QStringLiteral("正视 +Z"));
    addViewButton(ViewerCameraDirection::NegativeZ,
                  QStringLiteral(":/sonic/icons/nZ.png"),
                  QStringLiteral("正视 -Z"));
    addViewButton(ViewerCameraDirection::PositiveX,
                  QStringLiteral(":/sonic/icons/pX.png"),
                  QStringLiteral("正视 +X"));
    addViewButton(ViewerCameraDirection::NegativeX,
                  QStringLiteral(":/sonic/icons/nX.png"),
                  QStringLiteral("正视 -X"));
    addViewButton(ViewerCameraDirection::PositiveY,
                  QStringLiteral(":/sonic/icons/pY.png"),
                  QStringLiteral("正视 +Y"));
    addViewButton(ViewerCameraDirection::NegativeY,
                  QStringLiteral(":/sonic/icons/nY.png"),
                  QStringLiteral("正视 -Y"));
    for (QToolButton* button : findChildren<QToolButton*>()) {
        if (button->property("cameraDirection").isValid()) {
            viewGroup->addButton(button);
        }
    }

    auto* projectionGroup = new QButtonGroup(this);
    projectionGroup->setExclusive(true);
    perspectiveOnButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/perspective-on.png"),
        QStringLiteral("开启透视投影"), true);
    perspectiveOffButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/perspective-off.png"),
        QStringLiteral("关闭透视投影，使用正交投影"), true);
    projectionGroup->addButton(perspectiveOnButton_);
    projectionGroup->addButton(perspectiveOffButton_);
    secondRowLayout_->addWidget(perspectiveOnButton_);
    secondRowLayout_->addWidget(perspectiveOffButton_);
    registerVtkControl(perspectiveOnButton_);
    registerVtkControl(perspectiveOffButton_);
    registerMeshControl(perspectiveOnButton_);
    registerMeshControl(perspectiveOffButton_);

    resetButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/mesh.png"),
        QStringLiteral("适应窗口"));
    secondRowLayout_->addWidget(resetButton_);
    registerVtkControl(resetButton_);
    registerMeshControl(resetButton_);

    secondRowLayout_->addSpacing(4);
    playButton_ = makeToolButton(QStringLiteral(":/sonic/icons/play.png"),
                                 QStringLiteral("P"),
                                 QStringLiteral("连续播放时间步"), true);
    previousTimeButton_ = makeToolButton(
        QStringLiteral(":/sonic/icons/previous.png"),
        QStringLiteral("<"), QStringLiteral("上一时间步"));
    nextTimeButton_ = makeToolButton(QStringLiteral(":/sonic/icons/next.png"),
                                     QStringLiteral(">"),
                                     QStringLiteral("下一时间步"));
    firstTimeButton_ = makeToolButton(QStringLiteral(":/sonic/icons/first.png"),
                                      QStringLiteral("|<"),
                                      QStringLiteral("回到起始时间"));
    lastTimeButton_ = makeToolButton(QStringLiteral(":/sonic/icons/final.png"),
                                     QStringLiteral(">|"),
                                     QStringLiteral("跳到最后时间"));
    secondRowLayout_->addWidget(playButton_);
    secondRowLayout_->addWidget(previousTimeButton_);
    secondRowLayout_->addWidget(nextTimeButton_);
    secondRowLayout_->addWidget(firstTimeButton_);
    secondRowLayout_->addWidget(lastTimeButton_);

    playbackRateSlider_ = new QSlider(Qt::Horizontal, this);
    playbackRateSlider_->setObjectName(QStringLiteral("playbackRateSlider"));
    playbackRateSlider_->setRange(playbackSliderMinimum,
                                  playbackSliderMaximum);
    playbackRateSlider_->setValue(sliderPositionForPlaybackRate(24));
    playbackRateSlider_->setFixedWidth(112);
    playbackRateSlider_->setToolTip(
        QStringLiteral("动画帧率：指数调节 1–999 fps，常用的 1–60 fps 占据更多行程"));
    playbackRateLabel_ = new QLabel(QStringLiteral("24 fps"), this);
    playbackRateLabel_->setObjectName(QStringLiteral("playbackRateLabel"));
    playbackRateLabel_->setMinimumWidth(48);
    playbackRateLabel_->setAlignment(Qt::AlignCenter);
    secondRowLayout_->addWidget(playbackRateSlider_);
    secondRowLayout_->addWidget(playbackRateLabel_);
    registerVtkControl(playButton_);
    registerVtkControl(previousTimeButton_);
    registerVtkControl(nextTimeButton_);
    registerVtkControl(firstTimeButton_);
    registerVtkControl(lastTimeButton_);
    registerMeshControl(playButton_);
    registerMeshControl(previousTimeButton_);
    registerMeshControl(nextTimeButton_);
    registerMeshControl(firstTimeButton_);
    registerMeshControl(lastTimeButton_);
    registerVtkControl(playbackRateSlider_);
    registerVtkControl(playbackRateLabel_);
    registerMeshControl(playbackRateSlider_);
    registerMeshControl(playbackRateLabel_);

    timeLabel_ = new QLabel(QStringLiteral("t = -"), this);
    timeLabel_->setObjectName(QStringLiteral("viewerTimeLabel"));
    timeLabel_->setMinimumWidth(118);
    timeLabel_->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
    secondRowLayout_->addWidget(timeLabel_);
    registerVtkControl(timeLabel_);
    registerMeshControl(timeLabel_);

    auto* geometryToolGroup = new QButtonGroup(this);
    geometryToolGroup->setExclusive(true);
    addGeometryTool(ViewerGeometryTool::Vertex,
                    QStringLiteral(":/sonic/icons/vertex.png"),
                    QStringLiteral("顶点"));
    addGeometryTool(ViewerGeometryTool::Line,
                    QStringLiteral(":/sonic/icons/line.png"),
                    QStringLiteral("直线"));
    addGeometryTool(ViewerGeometryTool::CenterArc,
                    QStringLiteral(":/sonic/icons/arc.png"),
                    QStringLiteral("以中心画圆弧"));
    addGeometryTool(ViewerGeometryTool::LengthDimension,
                    QStringLiteral(":/sonic/icons/length.png"),
                    QStringLiteral("标识长度"));
    addGeometryTool(ViewerGeometryTool::AngleDimension,
                    QStringLiteral(":/sonic/icons/angle.png"),
                    QStringLiteral("标识角度"));
    addGeometryTool(ViewerGeometryTool::TrimSegment,
                    QStringLiteral(":/sonic/icons/trim.png"),
                    QStringLiteral("剪裁线段"));
    for (QToolButton* button : findChildren<QToolButton*>()) {
        if (button->property("geometryTool").isValid()) {
            geometryToolGroup->addButton(button);
        }
    }
    secondRowLayout_->addStretch();

    sliceButton_ = makeIconButton(
        QStringLiteral(":/sonic/icons/slice.png"),
        QStringLiteral("Slice：打开切片坐标面板"), true);
    thirdRowLayout_->addWidget(sliceButton_);
    registerVtkControl(sliceButton_);
    registerMeshControl(sliceButton_);

    extrudeButton_ = makeToolButton(QStringLiteral(":/sonic/icons/extrude.png"),
                                    QStringLiteral("E"),
                                    QStringLiteral("拉伸"));
    revolveButton_ = makeToolButton(QStringLiteral(":/sonic/icons/revolve.png"),
                                    QStringLiteral("R"),
                                    QStringLiteral("旋转扫描"));
    thirdRowLayout_->addWidget(extrudeButton_);
    thirdRowLayout_->addWidget(revolveButton_);
    registerGeometryControl(extrudeButton_);
    registerGeometryControl(revolveButton_);
    thirdRowLayout_->addStretch();

    connect(sliceButton_, &QToolButton::toggled,
            this, &ViewerToolbar::sliceToggled);
    connect(scalarSelector_, &QComboBox::currentTextChanged,
            this, &ViewerToolbar::scalarChanged);
    connect(surfaceButton_, &QToolButton::clicked, this, [this] {
        emit displayModeChanged(ViewerDisplayMode::PhysicalSurface);
    });
    connect(meshButton_, &QToolButton::clicked, this, [this] {
        emit displayModeChanged(ViewerDisplayMode::Mesh);
    });
    connect(perspectiveOnButton_, &QToolButton::clicked, this, [this] {
        emit perspectiveChanged(true);
    });
    connect(perspectiveOffButton_, &QToolButton::clicked, this, [this] {
        emit perspectiveChanged(false);
    });
    connect(resetButton_, &QToolButton::clicked,
            this, &ViewerToolbar::resetCameraRequested);
    connect(playButton_, &QToolButton::toggled,
            this, &ViewerToolbar::playbackToggled);
    connect(previousTimeButton_, &QToolButton::clicked,
            this, &ViewerToolbar::previousTimeStepRequested);
    connect(nextTimeButton_, &QToolButton::clicked,
            this, &ViewerToolbar::nextTimeStepRequested);
    connect(firstTimeButton_, &QToolButton::clicked,
            this, &ViewerToolbar::firstTimeStepRequested);
    connect(lastTimeButton_, &QToolButton::clicked,
            this, &ViewerToolbar::lastTimeStepRequested);
    connect(playbackRateSlider_, &QSlider::valueChanged, this,
            [this](int sliderPosition) {
        const int framesPerSecond =
            playbackRateForSliderPosition(sliderPosition);
        playbackRateLabel_->setText(
            QStringLiteral("%1 fps").arg(framesPerSecond));
        if (playbackRateHandler_) playbackRateHandler_(framesPerSecond);
    });
    connect(extrudeButton_, &QToolButton::clicked, this, [this] {
        bool accepted = false;
        const double length = QInputDialog::getDouble(
            this, QStringLiteral("拉伸"),
            QStringLiteral("拉伸长度"), 1.0, -1.0e9, 1.0e9, 4,
            &accepted);
        if (accepted) emit solidExtrudeRequested(length);
    });
    connect(revolveButton_, &QToolButton::clicked, this, [this] {
        bool accepted = false;
        const double angle = QInputDialog::getDouble(
            this, QStringLiteral("旋转扫描"),
            QStringLiteral("旋转角度"), 360.0, -3600.0, 3600.0, 2,
            &accepted);
        if (accepted) emit solidRevolveRequested(angle);
    });

    setDisplayMode(ViewerDisplayMode::PhysicalSurface);
    setPerspectiveEnabled(true);
    setViewDirection(ViewerCameraDirection::PositiveZ);
    setTimeFrameState(-1, 0, QString());
    setWorkspaceMode(ViewerWorkspaceMode::Geometry);
}

void ViewerToolbar::setVtkAvailable(bool available)
{
    for (QWidget* control : std::as_const(vtkControls_)) {
        if (control) control->setEnabled(available);
    }
}

void ViewerToolbar::setNavigationControls(QWidget* moduleSelector,
                                          QWidget* viewModeSwitch)
{
    if (!firstRowLayout_) return;
    if (moduleSelector) {
        moduleSelector->setParent(this);
        firstRowLayout_->insertWidget(0, moduleSelector);
    }
    if (viewModeSwitch) {
        viewModeSwitch->setParent(this);
        const int index = moduleSelector ? 1 : 0;
        firstRowLayout_->insertWidget(index, viewModeSwitch);
    }
}

void ViewerToolbar::setWorkspaceMode(ViewerWorkspaceMode mode)
{
    const bool meshMode = mode == ViewerWorkspaceMode::Mesh;
    for (QWidget* control : std::as_const(meshControls_)) {
        if (control) control->setVisible(meshMode);
    }
    for (QWidget* control : std::as_const(geometryControls_)) {
        if (control) control->setVisible(!meshMode);
    }
}

void ViewerToolbar::setScalarChoices(const QStringList& scalars,
                                     const QString& preferredScalar)
{
    const QSignalBlocker blocker(scalarSelector_);
    scalarSelector_->clear();
    scalarSelector_->addItems(scalars);

    int selected = -1;
    for (int i = 0; i < scalarSelector_->count(); ++i) {
        if (scalarSelector_->itemText(i).compare(preferredScalar,
                                                 Qt::CaseInsensitive) == 0) {
            selected = i;
            break;
        }
    }
    if (selected < 0 && !scalars.isEmpty()) selected = 0;
    if (selected >= 0) scalarSelector_->setCurrentIndex(selected);
    scalarSelector_->setEnabled(!scalars.isEmpty());
    update();
}

void ViewerToolbar::setSliceEnabled(bool enabled)
{
    const QSignalBlocker blocker(sliceButton_);
    sliceButton_->setChecked(enabled);
}

void ViewerToolbar::setDisplayMode(ViewerDisplayMode mode)
{
    const QSignalBlocker surfaceBlocker(surfaceButton_);
    const QSignalBlocker meshBlocker(meshButton_);
    surfaceButton_->setChecked(mode == ViewerDisplayMode::PhysicalSurface);
    meshButton_->setChecked(mode == ViewerDisplayMode::Mesh);
}

void ViewerToolbar::setPerspectiveEnabled(bool enabled)
{
    const QSignalBlocker onBlocker(perspectiveOnButton_);
    const QSignalBlocker offBlocker(perspectiveOffButton_);
    perspectiveOnButton_->setChecked(enabled);
    perspectiveOffButton_->setChecked(!enabled);
}

void ViewerToolbar::setViewDirection(ViewerCameraDirection direction)
{
    for (QToolButton* button : findChildren<QToolButton*>()) {
        const QVariant value = button->property("cameraDirection");
        if (!value.isValid()) continue;
        const auto buttonDirection =
            static_cast<ViewerCameraDirection>(value.toInt());
        const QSignalBlocker blocker(button);
        button->setChecked(buttonDirection == direction);
    }
}

void ViewerToolbar::setTimeFrameState(int currentIndex,
                                      int frameCount,
                                      const QString& physicalTime)
{
    const bool hasFrames = frameCount > 0 && currentIndex >= 0
        && currentIndex < frameCount;
    const bool multiFrame = frameCount > 1;
    const QString time = formatPhysicalTime(physicalTime);
    timeLabel_->setText(hasFrames
        ? QStringLiteral("t = %1  %2/%3")
              .arg(time)
              .arg(currentIndex + 1)
              .arg(frameCount)
        : QStringLiteral("t = -"));
    playButton_->setEnabled(multiFrame);
    previousTimeButton_->setEnabled(multiFrame);
    nextTimeButton_->setEnabled(multiFrame);
    firstTimeButton_->setEnabled(multiFrame);
    lastTimeButton_->setEnabled(multiFrame);
    playbackRateSlider_->setEnabled(multiFrame);
    playbackRateLabel_->setEnabled(multiFrame);
}

void ViewerToolbar::setPlaybackEnabled(bool enabled)
{
    const QSignalBlocker blocker(playButton_);
    playButton_->setChecked(enabled);
    playButton_->setIcon(QIcon(enabled
        ? QStringLiteral(":/sonic/icons/pause.png")
        : QStringLiteral(":/sonic/icons/play.png")));
}

void ViewerToolbar::setPlaybackRateHandler(
    std::function<void(int)> handler)
{
    playbackRateHandler_ = std::move(handler);
}

QToolButton* ViewerToolbar::makeIconButton(const QString& iconPath,
                                           const QString& tooltip,
                                           bool checkable)
{
    auto* button = new QToolButton(this);
    button->setIcon(QIcon(iconPath));
    button->setIconSize(QSize(16, 16));
    button->setFixedSize(28, 28);
    button->setToolTip(tooltip);
    button->setCheckable(checkable);
    button->setAutoRaise(false);
    button->setFocusPolicy(Qt::NoFocus);
    button->setAttribute(Qt::WA_StyledBackground, true);
    return button;
}

QToolButton* ViewerToolbar::makeToolButton(const QString& iconPath,
                                           const QString& fallbackText,
                                           const QString& tooltip,
                                           bool checkable)
{
    auto* button = makeIconButton(iconPath, tooltip, checkable);
    button->setText(fallbackText);
    button->setToolButtonStyle(QIcon(iconPath).isNull()
        ? Qt::ToolButtonTextOnly : Qt::ToolButtonIconOnly);
    return button;
}

void ViewerToolbar::addViewButton(ViewerCameraDirection direction,
                                  const QString& iconPath,
                                  const QString& tooltip)
{
    QToolButton* button = makeIconButton(iconPath, tooltip);
    button->setCheckable(true);
    button->setProperty("cameraDirection", static_cast<int>(direction));
    secondRowLayout_->addWidget(button);
    registerVtkControl(button);
    registerMeshControl(button);
    connect(button, &QToolButton::clicked, this, [this, direction] {
        emit viewDirectionRequested(direction);
    });
}

void ViewerToolbar::addGeometryTool(ViewerGeometryTool tool,
                                    const QString& iconPath,
                                    const QString& tooltip)
{
    QToolButton* button = makeToolButton(iconPath, QString(), tooltip, true);
    button->setProperty("geometryTool", static_cast<int>(tool));
    if (tool == ViewerGeometryTool::Vertex) button->setChecked(true);
    secondRowLayout_->addWidget(button);
    registerGeometryControl(button);
    connect(button, &QToolButton::clicked, this, [this, tool] {
        emit geometryToolChanged(tool);
    });
}

void ViewerToolbar::registerVtkControl(QWidget* control)
{
    if (!control || vtkControls_.contains(control)) return;
    vtkControls_.append(control);
}

void ViewerToolbar::registerMeshControl(QWidget* control)
{
    if (!control || meshControls_.contains(control)) return;
    meshControls_.append(control);
}

void ViewerToolbar::registerGeometryControl(QWidget* control)
{
    if (!control || geometryControls_.contains(control)) return;
    geometryControls_.append(control);
}

} // namespace SF::GUI
