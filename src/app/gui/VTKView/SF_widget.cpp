/// @file SF_widget.cpp
/// @brief GUI 内嵌 VTK 结果视图实现。

#include "SF_widget.h"

#include "SF_toolbar.h"

#include <QAbstractSpinBox>
#include <QApplication>
#include <QCursor>
#include <QDir>
#include <QDoubleSpinBox>
#include <QEvent>
#include <QFile>
#include <QFileSystemWatcher>
#include <QFileInfo>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QMap>
#include <QHash>
#include <QMouseEvent>
#include <QPushButton>
#include <QRegularExpression>
#include <QSignalBlocker>
#include <QStackedLayout>
#include <QSet>
#include <QTimer>
#include <QVBoxLayout>
#include <QVector3D>
#include <QXmlStreamReader>

#ifdef SONIC_GUI_WITH_VTK
#include "SF_interactorStyle.h"

#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkCutter.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGeometryFilter.h>
#include <vtkFieldData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkOpenFOAMReader.h>
#include <vtkNew.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkThreshold.h>
#include <vtkIntArray.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkTextProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif

#include <cmath>
#include <algorithm>
#include <utility>

/// @file SF_widget.cpp
/// @brief GUI 右上角 VTK 3D 预览窗口实现。
///
/// 本文件只负责可视化层：读取 VTK/STL 输出、提取 surface、按标量着色、
/// 切片、相机控制和结果目录自动刷新。求解器数据生成仍在 CLI/solver 侧完成。

namespace SF::GUI {

#ifdef SONIC_GUI_WITH_VTK
namespace {

/// @brief 复制任意 vtkDataSet，避免 reader 生命周期结束后输出失效。
vtkSmartPointer<vtkDataSet> copyDataSet(vtkDataSet* input)
{
    vtkSmartPointer<vtkDataSet> output;
    if (!input) return output;
    output.TakeReference(input->NewInstance());
    output->DeepCopy(input);
    return output;
}

/// @brief 复制 vtkPolyData，统一交给 mapper 使用。
vtkSmartPointer<vtkPolyData> copyPolyData(vtkPolyData* input)
{
    auto output = vtkSmartPointer<vtkPolyData>::New();
    if (input) output->DeepCopy(input);
    return output;
}

/// @brief 从结构/非结构网格中提取外表面。
///
/// VTK mapper 最终只接收 `vtkPolyData`，因此 `.vts/.vtu` 读入后先通过
/// `vtkDataSetSurfaceFilter` 转为可渲染 surface。
vtkSmartPointer<vtkPolyData> surfaceFromDataSet(vtkDataSet* input)
{
    vtkNew<vtkDataSetSurfaceFilter> surface;
    surface->SetInputData(input);
    surface->Update();
    return copyPolyData(surface->GetOutput());
}

/// @brief 把多块 VTM 合并成一个几何 surface。
vtkSmartPointer<vtkPolyData> surfaceFromAlgorithm(vtkAlgorithmOutput* outputPort)
{
    vtkNew<vtkCompositeDataGeometryFilter> geometry;
    geometry->SetInputConnection(outputPort);
    geometry->Update();
    return copyPolyData(geometry->GetOutput());
}

/// @brief 严格判断两个时间步是否使用完全相同的点坐标和单元连接。
///
/// 只有确认拓扑完全一致才复用已有 actor/mapper；移动网格或拓扑变化会显式
/// 返回 false 并走完整重建，不能为了播放速度错误沿用旧几何。
bool sameGeometry(vtkDataSet* previous, vtkDataSet* next)
{
    if (!previous || !next
        || QString::fromLatin1(previous->GetClassName())
               != QString::fromLatin1(next->GetClassName())
        || previous->GetNumberOfPoints() != next->GetNumberOfPoints()
        || previous->GetNumberOfCells() != next->GetNumberOfCells()) {
        return false;
    }

    double previousPoint[3];
    double nextPoint[3];
    for (vtkIdType point = 0; point < previous->GetNumberOfPoints(); ++point) {
        previous->GetPoint(point, previousPoint);
        next->GetPoint(point, nextPoint);
        if (previousPoint[0] != nextPoint[0]
            || previousPoint[1] != nextPoint[1]
            || previousPoint[2] != nextPoint[2]) {
            return false;
        }
    }

    vtkNew<vtkIdList> previousIds;
    vtkNew<vtkIdList> nextIds;
    for (vtkIdType cell = 0; cell < previous->GetNumberOfCells(); ++cell) {
        previous->GetCellPoints(cell, previousIds);
        next->GetCellPoints(cell, nextIds);
        if (previousIds->GetNumberOfIds() != nextIds->GetNumberOfIds()) {
            return false;
        }
        for (vtkIdType id = 0; id < previousIds->GetNumberOfIds(); ++id) {
            if (previousIds->GetId(id) != nextIds->GetId(id)) return false;
        }
    }
    return true;
}

/// @brief 保留目标网格的 points/cells，只替换本时间步物理数组。
void copyFrameArrays(vtkDataSet* target, vtkDataSet* source)
{
    target->GetPointData()->DeepCopy(source->GetPointData());
    target->GetCellData()->DeepCopy(source->GetCellData());
    target->GetFieldData()->DeepCopy(source->GetFieldData());
    target->Modified();
}

/// @brief 按 cell-data 整数 ID 提取可直接给 vtkPolyDataMapper 的高亮几何。
vtkSmartPointer<vtkPolyData> cellsWithId(vtkDataSet* input,
                                         const char* arrayName,
                                         int id)
{
    auto threshold = vtkSmartPointer<vtkThreshold>::New();
    threshold->SetInputData(input);
    threshold->SetInputArrayToProcess(0, 0, 0,
        vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayName);
    threshold->SetLowerThreshold(id);
    threshold->SetUpperThreshold(id);
    threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    threshold->Update();
    auto geometry = vtkSmartPointer<vtkGeometryFilter>::New();
    geometry->SetInputConnection(threshold->GetOutputPort());
    geometry->Update();
    return copyPolyData(geometry->GetOutput());
}

struct SfmPreview {
    vtkSmartPointer<vtkStructuredGrid> grid;
    QStringList setLabels;
};

QString vtkSetScalarName(const QString& setName)
{
    QString sanitized = setName.trimmed();
    sanitized.replace(QRegularExpression(QStringLiteral(R"([^A-Za-z0-9_]+)")),
                      QStringLiteral("_"));
    return QStringLiteral("set_%1").arg(
        sanitized.isEmpty() ? QStringLiteral("unnamed") : sanitized);
}

QString cleanSfmLine(QString line)
{
    const int comment = line.indexOf(QStringLiteral("//"));
    if (comment >= 0) line = line.left(comment);
    return line.trimmed();
}

bool readSfmFile(const QString& path,
                 int* nx,
                 int* ny,
                 int* nz,
                 QVector<QVector3D>* points,
                 QMap<QString, QVector<int>>* sets,
                 QString* error)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }

    QString section;
    bool readingPoints = false;
    while (!file.atEnd()) {
        QString line = cleanSfmLine(QString::fromUtf8(file.readLine()));
        if (line.isEmpty()) continue;
        if (line.startsWith(QLatin1Char('#'))) {
            section = line.mid(1).trimmed();
            readingPoints = section.compare(QStringLiteral("Point"),
                                            Qt::CaseInsensitive) == 0;
            if (sets
                && section.compare(QStringLiteral("Information"),
                                   Qt::CaseInsensitive) != 0
                && section.compare(QStringLiteral("Point"),
                                   Qt::CaseInsensitive) != 0
                && !section.isEmpty()) {
                (*sets)[section];
            }
            continue;
        }

        const QStringList tokens = line.split(
            QRegularExpression(QStringLiteral("\\s+")), Qt::SkipEmptyParts);
        if (section.compare(QStringLiteral("Information"),
                            Qt::CaseInsensitive) == 0) {
            if (tokens.size() >= 3) {
                if (nx) *nx = tokens.at(0).toInt();
                if (ny) *ny = tokens.at(1).toInt();
                if (nz) *nz = tokens.at(2).toInt();
            }
            continue;
        }
        if (readingPoints && points && tokens.size() >= 3) {
            points->push_back(QVector3D(
                static_cast<float>(tokens.at(0).toDouble()),
                static_cast<float>(tokens.at(1).toDouble()),
                static_cast<float>(tokens.at(2).toDouble())));
            continue;
        }
        if (sets
            && section.compare(QStringLiteral("Information"),
                               Qt::CaseInsensitive) != 0
            && section.compare(QStringLiteral("Point"),
                               Qt::CaseInsensitive) != 0
            && !section.isEmpty()) {
            QVector<int>& indices = (*sets)[section];
            for (const QString& token : tokens) {
                bool ok = false;
                const int index = token.toInt(&ok);
                if (ok) indices.push_back(index);
            }
        }
    }
    return true;
}

QStringList sfmSetCandidates(const QFileInfo& info)
{
    QStringList candidates;
    Q_UNUSED(info);
    return candidates;
}

SfmPreview readSfmPreview(const QString& path, QString* error)
{
    int nx = 0;
    int ny = 0;
    int nz = 0;
    QVector<QVector3D> points;
    QMap<QString, QVector<int>> sets;
    if (!readSfmFile(path, &nx, &ny, &nz, &points, &sets, error)) return {};

    for (const QString& candidate : sfmSetCandidates(QFileInfo(path))) {
        if (!QFileInfo::exists(candidate)
            || QFileInfo(candidate).absoluteFilePath()
                   == QFileInfo(path).absoluteFilePath()) {
            continue;
        }
        QMap<QString, QVector<int>> companionSets;
        int unusedNx = 0;
        int unusedNy = 0;
        int unusedNz = 0;
        QVector<QVector3D> unusedPoints;
        QString companionError;
        if (readSfmFile(candidate, &unusedNx, &unusedNy, &unusedNz,
                        &unusedPoints, &companionSets, &companionError)) {
            for (auto it = companionSets.cbegin();
                 it != companionSets.cend(); ++it) {
                sets[it.key()] = it.value();
            }
            break;
        }
    }

    if (nx <= 0 || ny <= 0 || nz <= 0 || points.isEmpty()) {
        if (error) *error = QStringLiteral("SFM 文件缺少有效网格点");
        return {};
    }
    const qsizetype expected = static_cast<qsizetype>(nx) * ny * nz;
    if (points.size() != expected) {
        if (error) {
            *error = QStringLiteral("SFM 点数与 #Information 不一致：%1 != %2")
                .arg(points.size())
                .arg(expected);
        }
        return {};
    }

    vtkNew<vtkPoints> vtkPoints;
    vtkPoints->SetNumberOfPoints(points.size());
    for (qsizetype i = 0; i < points.size(); ++i) {
        const QVector3D& point = points.at(i);
        vtkPoints->SetPoint(i, point.x(), point.y(), point.z());
    }

    auto grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(nx, ny, nz);
    grid->SetPoints(vtkPoints);

    if (!sets.isEmpty()) {
        QVector<QString> orderedNames;
        if (sets.contains(QStringLiteral("all"))) {
            orderedNames.push_back(QStringLiteral("all"));
        }
        for (auto it = sets.cbegin(); it != sets.cend(); ++it) {
            if (it.key().compare(QStringLiteral("all"),
                                 Qt::CaseInsensitive) != 0) {
                orderedNames.push_back(it.key());
            }
        }

        QVector<int> ids(points.size(), 0);
        QStringList labels;
        labels << QStringLiteral("0:none");
        for (int i = 0; i < orderedNames.size(); ++i) {
            const QString& name = orderedNames.at(i);
            const int id = i + 1;
            labels << QStringLiteral("%1:%2").arg(id).arg(name);
            for (int index : sets.value(name)) {
                if (index >= 0 && index < ids.size()) ids[index] = id;
            }
        }

        vtkNew<vtkIntArray> setArray;
        setArray->SetName("SetID");
        setArray->SetNumberOfComponents(1);
        setArray->SetNumberOfTuples(ids.size());
        for (qsizetype i = 0; i < ids.size(); ++i) {
            setArray->SetValue(i, ids.at(i));
        }
        grid->GetPointData()->AddArray(setArray);
        grid->GetPointData()->SetActiveScalars("SetID");

        for (const QString& name : orderedNames) {
            QVector<int> mask(points.size(), 0);
            for (int index : sets.value(name)) {
                if (index >= 0 && index < mask.size()) mask[index] = 1;
            }
            vtkNew<vtkIntArray> setMask;
            const QByteArray arrayName = vtkSetScalarName(name).toUtf8();
            setMask->SetName(arrayName.constData());
            setMask->SetNumberOfComponents(1);
            setMask->SetNumberOfTuples(mask.size());
            for (qsizetype i = 0; i < mask.size(); ++i) {
                setMask->SetValue(i, mask.at(i));
            }
            grid->GetPointData()->AddArray(setMask);
        }
        return {grid, labels};
    }
    return {grid, {}};
}

/// @brief 把切片轴枚举转成数组下标。
int axisIndex(ViewerSliceAxis axis)
{
    switch (axis) {
        case ViewerSliceAxis::Y: return 1;
        case ViewerSliceAxis::Z: return 2;
        case ViewerSliceAxis::X:
        default: return 0;
    }
}

/// @brief 把数组下标转回切片轴枚举。
ViewerSliceAxis axisFromIndex(int index)
{
    switch (index) {
        case 1: return ViewerSliceAxis::Y;
        case 2: return ViewerSliceAxis::Z;
        default: return ViewerSliceAxis::X;
    }
}

/// @brief 选择包围盒最长方向作为初始切片法向。
///
/// 这样第一次打开 Slice 时，切片面通常落在模型主要尺度方向的中间位置。
ViewerSliceAxis longestAxis(const double bounds[6])
{
    const double length[3] = {
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4]
    };
    int axis = 0;
    if (length[1] > length[axis]) axis = 1;
    if (length[2] > length[axis]) axis = 2;
    return axisFromIndex(axis);
}

/// @brief 计算 VTK bounds 的中心点。
void boundsCenter(const double bounds[6], double center[3])
{
    center[0] = 0.5 * (bounds[0] + bounds[1]);
    center[1] = 0.5 * (bounds[2] + bounds[3]);
    center[2] = 0.5 * (bounds[4] + bounds[5]);
}

/// @brief 生成沿指定坐标轴的单位法向。
void normalFromAxis(ViewerSliceAxis axis, double normal[3])
{
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    normal[axisIndex(axis)] = 1.0;
}

/// @brief 根据输入范围给切片坐标框一个合适的步长。
double spinStep(double minimum, double maximum)
{
    const double span = maximum - minimum;
    return span > 0.0 ? span / 100.0 : 1.0e-6;
}

/// @brief 用 vtkCutter 按点和法向生成切片 PolyData。
///
/// `origin` 是切片面经过的点，`normal` 是切片面法向。这里不改变原始数据，
/// 而是把 cutter 输出缓存到 `sliceData`，后续由显示管线决定显示整面还是切片。
vtkSmartPointer<vtkPolyData> sliceFromDataSet(vtkDataSet* input,
                                              const double origin[3],
                                              const double normal[3])
{
    auto output = vtkSmartPointer<vtkPolyData>::New();
    if (!input || input->GetNumberOfPoints() == 0) return output;

    double planeNormal[3] = {normal[0], normal[1], normal[2]};
    if (vtkMath::Normalize(planeNormal) <= 1.0e-12) {
        planeNormal[0] = 1.0;
        planeNormal[1] = 0.0;
        planeNormal[2] = 0.0;
    }

    vtkNew<vtkPlane> plane;
    plane->SetOrigin(origin);
    plane->SetNormal(planeNormal);

    vtkNew<vtkCutter> cutter;
    cutter->SetInputData(input);
    cutter->SetCutFunction(plane);
    cutter->Update();
    output->DeepCopy(cutter->GetOutput());
    return output;
}

/// @brief 收集可用于物理量着色的一维点标量数组名。
///
/// 目前工具栏只暴露单分量 PointData；`vtkGhostType` 是 VTK 内部数组，不显示给用户。
QStringList scalarNames(vtkPolyData* data)
{
    QStringList names;
    if (!data || !data->GetPointData()) return names;
    vtkPointData* pointData = data->GetPointData();
    for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
        vtkDataArray* array = pointData->GetArray(i);
        if (!array || array->GetNumberOfComponents() != 1) continue;
        const char* name = array->GetName();
        if (!name || !*name) continue;
        const QString field = QString::fromUtf8(name);
        if (field.compare(QStringLiteral("vtkGhostType"),
                          Qt::CaseInsensitive) == 0) {
            continue;
        }
        names.push_back(field);
    }
    names.removeDuplicates();
    return names;
}

/// @brief 按用户偏好选择标量，默认优先 Pressure。
QString chooseScalar(const QStringList& names, const QString& preferred)
{
    for (const QString& name : names) {
        if (name.compare(preferred, Qt::CaseInsensitive) == 0) return name;
    }
    for (const QString& name : names) {
        if (name.compare(QStringLiteral("Pressure"),
                         Qt::CaseInsensitive) == 0) {
            return name;
        }
    }
    return names.isEmpty() ? QString() : names.first();
}

/// @brief 在当前 PolyData 中按名称查找点标量数组。
vtkDataArray* findScalarArray(vtkPolyData* data, const QString& scalarName)
{
    if (!data || scalarName.isEmpty() || !data->GetPointData()) return nullptr;
    vtkPointData* pointData = data->GetPointData();
    for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
        vtkDataArray* array = pointData->GetArray(i);
        if (!array || !array->GetName()) continue;
        if (QString::fromUtf8(array->GetName()).compare(
                scalarName, Qt::CaseInsensitive) == 0) {
            return array;
        }
    }
    return nullptr;
}

/// @brief 从 actor bounds 得到相机对焦中心和包围半径。
///
/// 视角按钮、正交投影缩放都用这个半径保证模型完整落在视野内。
void boundsCenterAndRadius(vtkActor* actor, double center[3], double* radius)
{
    double bounds[6] = {-0.5, 0.5, -0.5, 0.5, -0.5, 0.5};
    if (actor) actor->GetBounds(bounds);
    center[0] = 0.5 * (bounds[0] + bounds[1]);
    center[1] = 0.5 * (bounds[2] + bounds[3]);
    center[2] = 0.5 * (bounds[4] + bounds[5]);
    const double dx = bounds[1] - bounds[0];
    const double dy = bounds[3] - bounds[2];
    const double dz = bounds[5] - bounds[4];
    *radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
    if (*radius <= 1.0e-12) *radius = 1.0;
}

} // namespace

class VTKView::Impl {
public:
    /// @brief Qt/VTK OpenGL 画布。
    QVTKOpenGLNativeWidget* widget = nullptr;
    /// @brief VTK 渲染窗口和 renderer。
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    vtkNew<vtkRenderer> renderer;
    /// @brief 当前显示对象的 mapper/actor。
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    /// @brief Geometry 工作区的持久 actor；编辑只更新 mapper 输入数据。
    vtkSmartPointer<vtkPolyDataMapper> geometryLineMapper;
    vtkSmartPointer<vtkActor> geometryLineActor;
    vtkSmartPointer<vtkPolyDataMapper> geometryPatchMapper;
    vtkSmartPointer<vtkActor> geometryPatchActor;
    vtkSmartPointer<vtkPolyDataMapper> geometryDefaultFaceMapper;
    vtkSmartPointer<vtkActor> geometryDefaultFaceActor;
    bool geometryDataReady = false;
    bool geometrySceneVisible = false;
    /// @brief 当前标量着色的 color bar。
    vtkNew<vtkScalarBarActor> scalarBar;
    /// @brief 原始数据、完整外表面、切片面和当前实际显示的数据。
    vtkSmartPointer<vtkDataSet> sourceData;
    vtkSmartPointer<vtkPolyData> surfaceData;
    vtkSmartPointer<vtkPolyData> sliceData;
    vtkSmartPointer<vtkPolyData> displayData;
    /// @brief 高亮 actor / mapper。
    vtkSmartPointer<vtkActor> highlightActor;
    vtkSmartPointer<vtkPolyDataMapper> highlightMapper;
    bool highlightActive = false;
    /// @brief geometry cell-data 数组（供 vtkThreshold 按 ID 筛选）。
    vtkSmartPointer<vtkPolyData> geometryLineData;
    vtkSmartPointer<vtkPolyData> geometryBlockData;
    vtkSmartPointer<vtkPolyData> geometryPatchData;
    /// @brief 物理量标量着色表。
    vtkNew<vtkLookupTable> lookupTable;
};
#endif

VTKView::VTKView(QWidget* parent)
    : VTKView(nullptr, parent)
{
}

VTKView::VTKView(ViewerToolbar* toolbar, QWidget* parent)
    : QWidget(parent)
{
    // QFileSystemWatcher 可能在一次求解输出中连续触发多次。
    // 用单次定时器合并事件，避免刚写到一半就抢读 PVD/VTK 文件。
    resultWatcher_ = new QFileSystemWatcher(this);
    resultReloadTimer_ = new QTimer(this);
    resultReloadTimer_->setSingleShot(true);
    resultReloadTimer_->setInterval(180);
    playbackTimer_ = new QTimer(this);
    playbackTimer_->setTimerType(Qt::PreciseTimer);
    playbackTimer_->setInterval(42);

    // 根布局 = VTK 画布区域。工具栏可由主窗口统一托管，也可内部创建。
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(0, 0, 0, 0);
    root->setSpacing(0);

    toolbar_ = toolbar ? toolbar : new ViewerToolbar(this);
    toolbar_->setVtkAvailable(vtkAvailable());
    toolbar_->setTimeFrameState(-1, 0, QString());
    toolbar_->setPerspectiveEnabled(true);
    toolbar_->setViewDirection(ViewerCameraDirection::PositiveZ);
    if (!toolbar) root->addWidget(toolbar_);

    viewerArea_ = new QWidget(this);
    viewerArea_->setObjectName(QStringLiteral("viewerArea"));
    auto* viewerStack = new QStackedLayout(viewerArea_);
    viewerStack->setContentsMargins(0, 0, 0, 0);
    // StackAll 让切片输入面板叠在 VTK 画布左上角，而不是切换页面。
    viewerStack->setStackingMode(QStackedLayout::StackAll);

#ifdef SONIC_GUI_WITH_VTK
    // VTK 渲染核心：一个 OpenGL widget、一个 render window、一个 renderer。
    impl_ = new Impl;
    impl_->widget = new QVTKOpenGLNativeWidget(viewerArea_);
    impl_->widget->setMouseTracking(false);
    impl_->widget->installEventFilter(this);
    impl_->window->AddRenderer(impl_->renderer);
    impl_->window->SetMultiSamples(0);
    impl_->renderer->SetPreserveColorBuffer(false);
    impl_->renderer->SetPreserveDepthBuffer(false);
    impl_->renderer->SetErase(true);
    impl_->widget->setRenderWindow(impl_->window);
    auto interactorStyle = createParaViewInteractorStyle();
    interactorStyle->SetDefaultRenderer(impl_->renderer);
    // 启动时就绑定当前 renderer，避免首次鼠标进入窗口时 VTK 状态不完整。
    interactorStyle->SetCurrentRenderer(impl_->renderer);
    impl_->widget->interactor()->SetInteractorStyle(interactorStyle);
    impl_->renderer->SetBackground(0.94, 0.91, 0.86);
    impl_->renderer->SetBackground2(0.99, 0.97, 0.93);
    impl_->renderer->GradientBackgroundOn();
    impl_->lookupTable->SetHueRange(0.66, 0.0);
    impl_->lookupTable->SetSaturationRange(0.86, 0.88);
    impl_->lookupTable->SetValueRange(0.92, 0.96);
    impl_->lookupTable->Build();
    impl_->scalarBar->SetLookupTable(impl_->lookupTable);
    impl_->scalarBar->SetNumberOfLabels(5);
    impl_->scalarBar->SetMaximumWidthInPixels(84);
    impl_->scalarBar->SetPosition(0.88, 0.16);
    impl_->scalarBar->SetWidth(0.08);
    impl_->scalarBar->SetHeight(0.68);
    impl_->scalarBar->GetTitleTextProperty()->SetColor(0.22, 0.20, 0.16);
    impl_->scalarBar->GetLabelTextProperty()->SetColor(0.22, 0.20, 0.16);
    impl_->scalarBar->GetTitleTextProperty()->SetFontSize(12);
    impl_->scalarBar->GetLabelTextProperty()->SetFontSize(10);
    viewerStack->addWidget(impl_->widget);

    // Slice 按钮打开的小坐标条。坐标框只负责输入切片平面经过的点，
    // 法向在点击“确定”时取当前相机朝向，贴近 ParaView 当前视角切片的使用习惯。
    slicePanel_ = new QWidget(viewerArea_);
    slicePanel_->setObjectName(QStringLiteral("slicePlanePanel"));
    slicePanel_->setFixedWidth(228);
    auto* sliceLayout = new QGridLayout(slicePanel_);
    sliceLayout->setContentsMargins(6, 5, 6, 5);
    sliceLayout->setHorizontalSpacing(4);
    sliceLayout->setVerticalSpacing(4);

    auto makeCoordinateInput = [this](const QString& name) {
        auto* input = new QDoubleSpinBox(slicePanel_);
        input->setObjectName(name);
        input->setDecimals(4);
        input->setButtonSymbols(QAbstractSpinBox::NoButtons);
        input->setKeyboardTracking(false);
        input->setFixedSize(54, 24);
        return input;
    };
    auto* xLabel = new QLabel(QStringLiteral("X"), slicePanel_);
    auto* yLabel = new QLabel(QStringLiteral("Y"), slicePanel_);
    auto* zLabel = new QLabel(QStringLiteral("Z"), slicePanel_);
    xLabel->setFixedWidth(10);
    yLabel->setFixedWidth(10);
    zLabel->setFixedWidth(10);
    sliceXInput_ = makeCoordinateInput(QStringLiteral("sliceXInput"));
    sliceYInput_ = makeCoordinateInput(QStringLiteral("sliceYInput"));
    sliceZInput_ = makeCoordinateInput(QStringLiteral("sliceZInput"));
    auto* applyButton = new QPushButton(QStringLiteral("确定"), slicePanel_);
    auto* cancelButton = new QPushButton(QStringLiteral("取消"), slicePanel_);
    applyButton->setObjectName(QStringLiteral("primaryButton"));
    applyButton->setFixedSize(42, 24);
    cancelButton->setFixedSize(42, 24);
    sliceLayout->addWidget(xLabel, 0, 0);
    sliceLayout->addWidget(sliceXInput_, 0, 1);
    sliceLayout->addWidget(yLabel, 0, 2);
    sliceLayout->addWidget(sliceYInput_, 0, 3);
    sliceLayout->addWidget(zLabel, 0, 4);
    sliceLayout->addWidget(sliceZInput_, 0, 5);
    sliceLayout->addWidget(cancelButton, 1, 3, 1, 2);
    sliceLayout->addWidget(applyButton, 1, 5);
    slicePanel_->hide();
    viewerStack->addWidget(slicePanel_);
    viewerStack->setAlignment(slicePanel_, Qt::AlignLeft | Qt::AlignTop);

    connect(applyButton, &QPushButton::clicked,
            this, &VTKView::applySliceFromPanel);
    connect(cancelButton, &QPushButton::clicked,
            this, &VTKView::closeSlicePanel);
    // 全局过滤空鼠标移动，避免普通划过 3D 窗口时触发相机动作。
    qApp->installEventFilter(this);
#else
    placeholder_ = new QLabel(
        QStringLiteral("VTK 预览尚未启用\n\n安装带 Qt 支持的 VTK 后重新配置项目，"
                       "这里会直接显示 .vts / .vtm 网格。"),
        viewerArea_);
    placeholder_->setAlignment(Qt::AlignCenter);
    placeholder_->setObjectName(QStringLiteral("viewerPlaceholder"));
    viewerStack->addWidget(placeholder_);
#endif
    root->addWidget(viewerArea_, 1);

    // 工具栏只发信号，真正的 VTK 状态切换都集中在 VTKView 内部完成。
    connect(toolbar_, &ViewerToolbar::resetCameraRequested,
            this, &VTKView::resetCamera);
    connect(toolbar_, &ViewerToolbar::sliceToggled,
            this, &VTKView::setSliceEnabled);
    connect(toolbar_, &ViewerToolbar::scalarChanged,
            this, &VTKView::setScalarField);
    connect(toolbar_, &ViewerToolbar::displayModeChanged,
            this, &VTKView::setDisplayMode);
    connect(toolbar_, &ViewerToolbar::perspectiveChanged,
            this, &VTKView::setPerspectiveEnabled);
    connect(toolbar_, &ViewerToolbar::viewDirectionRequested,
            this, &VTKView::setCameraDirection);
    connect(toolbar_, &ViewerToolbar::playbackToggled,
            this, &VTKView::setPlaybackEnabled);
    toolbar_->setPlaybackRateHandler(
        [this](int framesPerSecond) { setPlaybackRate(framesPerSecond); });
    connect(toolbar_, &ViewerToolbar::previousTimeStepRequested,
            this, &VTKView::showPreviousTimeStep);
    connect(toolbar_, &ViewerToolbar::nextTimeStepRequested,
            this, &VTKView::showNextTimeStep);
    connect(toolbar_, &ViewerToolbar::firstTimeStepRequested,
            this, &VTKView::showFirstTimeStep);
    connect(toolbar_, &ViewerToolbar::lastTimeStepRequested,
            this, &VTKView::showLastTimeStep);
    connect(resultWatcher_, &QFileSystemWatcher::directoryChanged,
            this, [this] { scheduleResultRefresh(); });
    connect(resultWatcher_, &QFileSystemWatcher::fileChanged,
            this, [this] { scheduleResultRefresh(); });
    connect(resultReloadTimer_, &QTimer::timeout,
            this, &VTKView::refreshResult);
    connect(playbackTimer_, &QTimer::timeout,
            this, &VTKView::advancePlayback);
}

VTKView::~VTKView()
{
#ifdef SONIC_GUI_WITH_VTK
    qApp->removeEventFilter(this);
    delete impl_;
#endif
}

bool VTKView::vtkAvailable() const
{
#ifdef SONIC_GUI_WITH_VTK
    return true;
#else
    return false;
#endif
}

bool VTKView::eventFilter(QObject* watched, QEvent* event)
{
#ifdef SONIC_GUI_WITH_VTK
    if (impl_ && impl_->widget && impl_->widget->isVisible()) {
        const bool cursorInViewer = impl_->widget->rect().contains(
            impl_->widget->mapFromGlobal(QCursor::pos()));
        // QVTK 在不同平台上会把 hover/tablet/native gesture 转成相机事件。
        // 这里吞掉“没有按键”的移动类事件，保留按住左/右键时的正常交互。
        if (event->type() == QEvent::MouseMove
            || event->type() == QEvent::HoverMove
            || event->type() == QEvent::TabletMove
            || event->type() == QEvent::NativeGesture) {
            if (!cursorInViewer) return QWidget::eventFilter(watched, event);
            if (event->type() == QEvent::MouseMove) {
                auto* mouseEvent = static_cast<QMouseEvent*>(event);
                if (mouseEvent->buttons() != Qt::NoButton) {
                    return QWidget::eventFilter(watched, event);
                }
            } else if (QApplication::mouseButtons() != Qt::NoButton) {
                return QWidget::eventFilter(watched, event);
            }
            return true;
        }
    }
#endif
    return QWidget::eventFilter(watched, event);
}

bool VTKView::loadBlockMeshGeometry(const QString& path, QString* error)
{
#ifdef SONIC_GUI_WITH_VTK
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }
    QString text = QString::fromUtf8(file.readAll());
    text.remove(QRegularExpression(QStringLiteral(R"(/\*[\s\S]*?\*/)")));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\n]*)")));

    auto listBody = [&text](const QString& key) {
        const QRegularExpression keyPattern(
            QStringLiteral(R"(\b%1\b\s*\()").arg(QRegularExpression::escape(key)));
        const QRegularExpressionMatch match = keyPattern.match(text);
        if (!match.hasMatch()) return QString();
        const int open = match.capturedEnd() - 1;
        int depth = 0;
        for (int i = open; i < text.size(); ++i) {
            if (text.at(i) == QLatin1Char('(')) ++depth;
            else if (text.at(i) == QLatin1Char(')') && --depth == 0)
                return text.mid(open + 1, i - open - 1);
        }
        return QString();
    };
    const QString number = QStringLiteral(
        R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)");
    const QString verticesBody = listBody(QStringLiteral("vertices"));
    const QRegularExpression vertex(QStringLiteral(
        R"(\(\s*(%1)\s+(%1)\s+(%1)\s*\))").arg(number));
    QVector<QVector3D> vertices;
    for (auto it = vertex.globalMatch(verticesBody); it.hasNext();) {
        const QRegularExpressionMatch match = it.next();
        vertices.push_back(QVector3D(match.captured(1).toFloat(),
                                     match.captured(2).toFloat(),
                                     match.captured(3).toFloat()));
    }
    if (vertices.isEmpty()) {
        if (error) *error = QStringLiteral("blockMeshDict 未包含有效 vertices");
        return false;
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (const QVector3D& vertexPoint : vertices) {
        points->InsertNextPoint(vertexPoint.x(), vertexPoint.y(), vertexPoint.z());
    }
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkCellArray>::New();
    auto blockFaces = vtkSmartPointer<vtkCellArray>::New();
    auto cellColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    cellColors->SetName("GeometryColor");
    cellColors->SetNumberOfComponents(3);
    const auto addLine = [&lines, &cellColors](vtkIdType first, vtkIdType second,
                                               unsigned char r, unsigned char g,
                                               unsigned char b) {
        vtkIdType ids[2] = {first, second};
        lines->InsertNextCell(2, ids);
        const unsigned char color[3] = {r, g, b};
        cellColors->InsertNextTypedTuple(color);
    };
    const auto indices = [](const QString& value) {
        QVector<int> result;
        const QRegularExpression integer(QStringLiteral(R"([-+]?\d+)"));
        for (auto it = integer.globalMatch(value); it.hasNext();)
            result.push_back(it.next().captured().toInt());
        return result;
    };
    const auto valid = [&vertices](int index) {
        return index >= 0 && index < vertices.size();
    };
    const auto canonicalFaceKey = [](QVector<int> faceVertices) {
        std::sort(faceVertices.begin(), faceVertices.end());
        QStringList values;
        values.reserve(faceVertices.size());
        for (int index : faceVertices) values.push_back(QString::number(index));
        return values.join(QLatin1Char(':'));
    };
    QHash<QString, int> blockFaceOccurrences;
    QHash<QString, QVector<int>> externalBlockFaces;
    QSet<QString> assignedPatchFaces;

    const QString blocksBody = listBody(QStringLiteral("blocks"));
    const QRegularExpression block(QStringLiteral(
        R"(\b(?:hex|wedge|prism)\s*\(([^)]*)\))"),
        QRegularExpression::CaseInsensitiveOption);
    static const int hexEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
        {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    static const int hexFaces[6][4] = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
        {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    auto blockIdArray = vtkSmartPointer<vtkIntArray>::New();
    blockIdArray->SetName("blockID");
    auto blockFaceIdArray = vtkSmartPointer<vtkIntArray>::New();
    blockFaceIdArray->SetName("blockID");
    int blockIdx = 0;
    for (auto it = block.globalMatch(blocksBody); it.hasNext(); ++blockIdx) {
        const QVector<int> blockVertices = indices(it.next().captured(1));
        if (blockVertices.size() != 8) continue;
        for (const auto& edge : hexEdges) {
            if (valid(blockVertices[edge[0]]) && valid(blockVertices[edge[1]])) {
                addLine(blockVertices[edge[0]], blockVertices[edge[1]], 180, 180, 180);
                blockIdArray->InsertNextValue(blockIdx);
            }
        }
        for (const auto& face : hexFaces) {
            QVector<int> faceVertices;
            faceVertices.reserve(4);
            for (int localIndex : face) faceVertices.push_back(blockVertices[localIndex]);
            if (std::any_of(faceVertices.cbegin(), faceVertices.cend(),
                            [&valid](int index) { return !valid(index); })) continue;
            const QString key = canonicalFaceKey(faceVertices);
            ++blockFaceOccurrences[key];
            externalBlockFaces.insert(key, faceVertices);
            vtkIdType ids[4] = {faceVertices[0], faceVertices[1],
                                faceVertices[2], faceVertices[3]};
            blockFaces->InsertNextCell(4, ids);
            blockFaceIdArray->InsertNextValue(blockIdx);
        }
    }

    const QString edgesBody = listBody(QStringLiteral("edges"));
    const QRegularExpression arc(QStringLiteral(
        R"(\barc\s+(\d+)\s+(\d+)\s*\(\s*(%1)\s+(%1)\s+(%1)\s*\))")
            .arg(number), QRegularExpression::CaseInsensitiveOption);
    for (auto it = arc.globalMatch(edgesBody); it.hasNext();) {
        const QRegularExpressionMatch match = it.next();
        const int first = match.captured(1).toInt();
        const int last = match.captured(2).toInt();
        if (!valid(first) || !valid(last)) continue;
        const QVector3D middle(match.captured(3).toFloat(), match.captured(4).toFloat(),
                               match.captured(5).toFloat());
        const QVector3D control = 2.0F * middle - 0.5F * (vertices[first] + vertices[last]);
        vtkIdType previous = first;
        for (int sample = 1; sample <= 32; ++sample) {
            const float t = static_cast<float>(sample) / 32.0F;
            const QVector3D point = (1.0F - t) * (1.0F - t) * vertices[first]
                + 2.0F * (1.0F - t) * t * control + t * t * vertices[last];
            const vtkIdType current = sample == 32
                ? static_cast<vtkIdType>(last)
                : points->InsertNextPoint(point.x(), point.y(), point.z());
            addLine(previous, current, 37, 99, 235);
            previous = current;
        }
    }

    auto patchIdArray = vtkSmartPointer<vtkIntArray>::New();
    patchIdArray->SetName("patchID");
    auto faceIdArray = vtkSmartPointer<vtkIntArray>::New();
    faceIdArray->SetName("faceID");
    const QString boundaryBody = listBody(QStringLiteral("boundary"));
    const QRegularExpression patchStart(
        QStringLiteral(R"(\b[A-Za-z_][A-Za-z0-9_.-]*\s*\{)"));
    int patchIdx = 0;
    int from = 0;
    while (true) {
        const QRegularExpressionMatch patch = patchStart.match(boundaryBody, from);
        if (!patch.hasMatch()) break;
        const int open = patch.capturedEnd() - 1;
        int depth = 0;
        int close = -1;
        for (int i = open; i < boundaryBody.size(); ++i) {
            if (boundaryBody.at(i) == QLatin1Char('{')) ++depth;
            else if (boundaryBody.at(i) == QLatin1Char('}') && --depth == 0) {
                close = i;
                break;
            }
        }
        if (close < 0) break;
        const QString patchText = boundaryBody.mid(open + 1, close - open - 1);
        int faceIdx = 0;
        const QRegularExpression facesKey(QStringLiteral(R"(\bfaces\b\s*\()"));
        const QRegularExpressionMatch facesMatch = facesKey.match(patchText);
        if (facesMatch.hasMatch()) {
            const int facesOpen = facesMatch.capturedEnd() - 1;
            int facesDepth = 0;
            int facesClose = -1;
            for (int i = facesOpen; i < patchText.size(); ++i) {
                if (patchText.at(i) == QLatin1Char('(')) ++facesDepth;
                else if (patchText.at(i) == QLatin1Char(')') && --facesDepth == 0) {
                    facesClose = i;
                    break;
                }
            }
            if (facesClose >= 0) {
                const QString faceText = patchText.mid(facesOpen + 1,
                                                       facesClose - facesOpen - 1);
                const QRegularExpression face(QStringLiteral(R"(\(([^()]*)\))"));
                for (auto fit = face.globalMatch(faceText); fit.hasNext();) {
                    const QVector<int> faceVertices = indices(fit.next().captured(1));
                    if (faceVertices.size() < 3
                        || std::any_of(faceVertices.cbegin(), faceVertices.cend(),
                                       [&valid](int index) { return !valid(index); })) continue;
                    QVector<vtkIdType> ids;
                    ids.reserve(faceVertices.size());
                    for (int index : faceVertices) ids.push_back(index);
                    faces->InsertNextCell(ids.size(), ids.data());
                    patchIdArray->InsertNextValue(patchIdx);
                    faceIdArray->InsertNextValue(faceIdx++);
                    assignedPatchFaces.insert(canonicalFaceKey(faceVertices));
                }
            }
        }
        ++patchIdx;
        from = close + 1;
    }

    polyData->SetPoints(points);
    polyData->SetLines(lines);
    polyData->GetCellData()->AddArray(blockIdArray);
    polyData->GetCellData()->SetScalars(cellColors);
    impl_->geometryLineData = polyData;
    auto blockData = vtkSmartPointer<vtkPolyData>::New();
    blockData->SetPoints(points);
    blockData->SetPolys(blockFaces);
    blockData->GetCellData()->AddArray(blockFaceIdArray);
    impl_->geometryBlockData = blockData;
    auto patchData = vtkSmartPointer<vtkPolyData>::New();
    patchData->SetPoints(points);
    patchData->GetCellData()->AddArray(patchIdArray);
    patchData->GetCellData()->AddArray(faceIdArray);
    patchData->SetPolys(faces);
    impl_->geometryPatchData = patchData;
    auto defaultFaces = vtkSmartPointer<vtkCellArray>::New();
    for (auto it = externalBlockFaces.cbegin(); it != externalBlockFaces.cend(); ++it) {
        if (blockFaceOccurrences.value(it.key()) != 1 || assignedPatchFaces.contains(it.key())) {
            continue;
        }
        QVector<vtkIdType> ids;
        ids.reserve(it.value().size());
        for (int index : it.value()) ids.push_back(index);
        defaultFaces->InsertNextCell(ids.size(), ids.data());
    }
    auto defaultFaceData = vtkSmartPointer<vtkPolyData>::New();
    defaultFaceData->SetPoints(points);
    defaultFaceData->SetPolys(defaultFaces);
    clearHighlight();
    const bool firstGeometryLoad = !impl_->geometryDataReady;
    if (!impl_->geometryLineMapper) {
        impl_->geometryLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        impl_->geometryLineMapper->SetScalarModeToUseCellData();
        impl_->geometryLineMapper->ScalarVisibilityOff();
        impl_->geometryLineActor = vtkSmartPointer<vtkActor>::New();
        impl_->geometryLineActor->SetMapper(impl_->geometryLineMapper);
        impl_->geometryLineActor->GetProperty()->SetColor(0.05, 0.05, 0.05);
        impl_->geometryLineActor->GetProperty()->SetOpacity(1.0);
        impl_->geometryLineActor->GetProperty()->SetLineWidth(1.5);
    }
    if (!impl_->geometryPatchMapper) {
        impl_->geometryPatchMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        impl_->geometryPatchMapper->ScalarVisibilityOff();
        impl_->geometryPatchActor = vtkSmartPointer<vtkActor>::New();
        impl_->geometryPatchActor->SetMapper(impl_->geometryPatchMapper);
        impl_->geometryPatchActor->GetProperty()->SetOpacity(1.0);
        impl_->geometryPatchActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    }
    if (!impl_->geometryDefaultFaceMapper) {
        impl_->geometryDefaultFaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        impl_->geometryDefaultFaceMapper->ScalarVisibilityOff();
        impl_->geometryDefaultFaceActor = vtkSmartPointer<vtkActor>::New();
        impl_->geometryDefaultFaceActor->SetMapper(impl_->geometryDefaultFaceMapper);
        impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(1.0);
        impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    }
    impl_->geometryLineMapper->SetInputData(polyData);
    impl_->geometryPatchMapper->SetInputData(patchData);
    impl_->geometryDefaultFaceMapper->SetInputData(defaultFaceData);
    impl_->geometryDataReady = true;
    impl_->geometryPatchActor->GetProperty()->SetOpacity(1.0);
    // 几何预览的面不依赖光源方向，否则无显式法向的背面会被误绘成黑色。
    impl_->geometryPatchActor->GetProperty()->SetAmbient(1.0);
    impl_->geometryPatchActor->GetProperty()->SetDiffuse(0.0);
    auto backface = vtkSmartPointer<vtkProperty>::New();
    // 面朝镜头时为橙色；背对镜头时为紫色，便于确认 patch 的朝向。
    impl_->geometryPatchActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    backface->SetColor(1.0, 1.0, 1.0);
    backface->SetOpacity(1.0);
    backface->SetAmbient(1.0);
    backface->SetDiffuse(0.0);
    impl_->geometryPatchActor->SetBackfaceProperty(backface);
    impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(1.0);
    impl_->geometryDefaultFaceActor->GetProperty()->SetAmbient(1.0);
    impl_->geometryDefaultFaceActor->GetProperty()->SetDiffuse(0.0);
    auto defaultBackface = vtkSmartPointer<vtkProperty>::New();
    defaultBackface->SetColor(1.0, 1.0, 1.0);
    defaultBackface->SetOpacity(1.0);
    defaultBackface->SetAmbient(1.0);
    defaultBackface->SetDiffuse(0.0);
    impl_->geometryDefaultFaceActor->SetBackfaceProperty(defaultBackface);
    sliceEnabled_ = false;
    if (impl_->geometrySceneVisible || !impl_->actor) {
        showGeometryScene();
        if (firstGeometryLoad) impl_->renderer->ResetCamera();
        impl_->renderer->ResetCameraClippingRange();
        impl_->window->Render();
    }
    emit statusMessage(QStringLiteral("Geometry 预览：%1 顶点").arg(vertices.size()));
    return true;
#else
    Q_UNUSED(path);
    if (error) *error = QStringLiteral("当前构建未启用 VTK，无法显示三维 Geometry 预览");
    return false;
#endif
}

bool VTKView::loadFile(const QString& path, QString* error)
{
    // PVD 是时间序列索引文件，真正的数据文件在 DataSet::file 属性里。
    const QFileInfo requested(path);
    if (requested.suffix().compare(QStringLiteral("pvd"),
                                   Qt::CaseInsensitive) == 0) {
        preferredPvd_ = requested.absoluteFilePath();
        return loadPvd(preferredPvd_, error);
    }

    setPlaybackEnabled(false);
    timeFrames_.clear();
    currentFrameIndex_ = -1;
    loadedPvdPath_.clear();
    preferredPvd_.clear();
    updateTimeControls();
    const bool loaded =
        loadDataFile(requested.absoluteFilePath(), error, true, false);
    if (loaded && visualizationLoadedHandler_) {
        visualizationLoadedHandler_(requested.absoluteFilePath());
    }
    return loaded;
}

bool VTKView::loadDataFile(const QString& path,
                           QString* error,
                           bool resetView,
                           bool fromTimeSeries)
{
#ifdef SONIC_GUI_WITH_VTK
    // 时间序列换帧必须严格保留用户相机；不能依赖 mapper/actor 更换过程中
    // VTK 的隐式状态，因为它可能重新计算距离并表现成画面逐帧缩小。
    double savedPosition[3] = {0.0, 0.0, 1.0};
    double savedFocalPoint[3] = {0.0, 0.0, 0.0};
    double savedViewUp[3] = {0.0, 1.0, 0.0};
    double savedClippingRange[2] = {0.01, 1000.0};
    double savedParallelScale = 1.0;
    double savedViewAngle = 30.0;
    int savedParallelProjection = 0;
    vtkCamera* savedCamera = impl_->renderer->GetActiveCamera();
    if (!resetView && savedCamera) {
        savedCamera->GetPosition(savedPosition);
        savedCamera->GetFocalPoint(savedFocalPoint);
        savedCamera->GetViewUp(savedViewUp);
        savedCamera->GetClippingRange(savedClippingRange);
        savedParallelScale = savedCamera->GetParallelScale();
        savedViewAngle = savedCamera->GetViewAngle();
        savedParallelProjection = savedCamera->GetParallelProjection();
    }
    const QFileInfo info(path);
    if (!info.isFile()) {
        if (error) *error = QStringLiteral("结果文件不存在：%1").arg(path);
        return false;
    }

    const QString suffix = info.suffix().toLower();
    vtkSmartPointer<vtkDataSet> sourceData;
    vtkSmartPointer<vtkPolyData> surfaceData;
    QStringList sfmSetLabels;

    // 不同 VTK 文件格式读进来的数据类型不同，这里统一转换到：
    // sourceData = 可切片的原始数据；surfaceData = 当前完整外表面显示数据。
    if (suffix == QStringLiteral("sfm")) {
        const SfmPreview preview = readSfmPreview(
            info.absoluteFilePath(), error);
        if (!preview.grid) return false;
        sourceData = preview.grid;
        surfaceData = surfaceFromDataSet(sourceData);
        selectedScalar_ = QStringLiteral("SetID");
        sfmSetLabels = preview.setLabels;
    } else if (suffix == QStringLiteral("vts")) {
        vtkNew<vtkXMLStructuredGridReader> reader;
        reader->SetFileName(info.absoluteFilePath().toUtf8().constData());
        reader->Update();
        sourceData = copyDataSet(reader->GetOutput());
        surfaceData = surfaceFromDataSet(sourceData);
    } else if (suffix == QStringLiteral("vtu")) {
        vtkNew<vtkXMLUnstructuredGridReader> reader;
        reader->SetFileName(info.absoluteFilePath().toUtf8().constData());
        reader->Update();
        sourceData = copyDataSet(reader->GetOutput());
        surfaceData = surfaceFromDataSet(sourceData);
    } else if (suffix == QStringLiteral("vtm")) {
        vtkNew<vtkXMLMultiBlockDataReader> reader;
        reader->SetFileName(info.absoluteFilePath().toUtf8().constData());
        reader->Update();
        surfaceData = surfaceFromAlgorithm(reader->GetOutputPort());
        sourceData = surfaceData;
    } else if (suffix == QStringLiteral("foam")) {
        vtkNew<vtkOpenFOAMReader> reader;
        reader->SetFileName(info.absoluteFilePath().toUtf8().constData());
        reader->SetCreateCellToPoint(true);
        reader->Update();
        surfaceData = surfaceFromAlgorithm(reader->GetOutputPort());
        sourceData = surfaceData;
    } else if (suffix == QStringLiteral("stl")) {
        vtkNew<vtkSTLReader> reader;
        reader->SetFileName(info.absoluteFilePath().toUtf8().constData());
        reader->Update();
        surfaceData = copyPolyData(reader->GetOutput());
        sourceData = surfaceData;
    } else {
        if (error) {
            *error = QStringLiteral(
                "仅支持 .sfm、.foam、.vts、.vtu、.vtm 和 .stl 文件");
        }
        return false;
    }

    const bool reuseStaticGeometry =
        fromTimeSeries
        && topologyPvdPath_ == loadedPvdPath_
        && impl_->actor && !impl_->geometrySceneVisible
        && sameGeometry(impl_->sourceData, sourceData)
        && sameGeometry(impl_->surfaceData, surfaceData);
    if (reuseStaticGeometry) {
        copyFrameArrays(impl_->sourceData, sourceData);
        copyFrameArrays(impl_->surfaceData, surfaceData);
        if (sliceEnabled_) updateSliceData();
        updateScalarChoices();
        updateVisualizationPipeline();
        loadedDataPath_ = info.absoluteFilePath();
        loadedDataModified_ = info.lastModified().toMSecsSinceEpoch();
        loadedDataSize_ = info.size();
        return true;
    }

    // Geometry 与 Mesh/Result 场景互斥。旧实现把几何面在每个时间步再次加入
    // renderer，与结果面共面叠加，旋转时会出现大面积 z-fighting 花屏。
    if (impl_->geometrySceneVisible) {
        if (impl_->highlightActor) impl_->renderer->RemoveActor(impl_->highlightActor);
        impl_->highlightActor = nullptr;
        impl_->highlightMapper = nullptr;
        if (impl_->geometryLineActor)
            impl_->renderer->RemoveActor(impl_->geometryLineActor);
        if (impl_->geometryPatchActor)
            impl_->renderer->RemoveActor(impl_->geometryPatchActor);
        if (impl_->geometryDefaultFaceActor)
            impl_->renderer->RemoveActor(impl_->geometryDefaultFaceActor);
        impl_->geometrySceneVisible = false;
        impl_->actor = nullptr;
        impl_->mapper = nullptr;
    }

    impl_->sourceData = sourceData;
    impl_->surfaceData = surfaceData;
    // 读取新数据后重置切片初值：原点取包围盒中心，法向取最长轴方向。
    impl_->sourceData->GetBounds(sliceBounds_);
    boundsCenter(sliceBounds_, sliceOrigin_);
    normalFromAxis(longestAxis(sliceBounds_), sliceNormal_);
    updateSliceData();
    // 先移除上一帧 actor，再替换智能指针。旧实现先覆盖指针，实际移除的是
    // 尚未加入 renderer 的新 actor，导致播放时 actor 持续累积。
    if (impl_->actor) impl_->renderer->RemoveActor(impl_->actor);
    impl_->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    impl_->actor = vtkSmartPointer<vtkActor>::New();
    impl_->actor->SetMapper(impl_->mapper);
    impl_->actor->GetProperty()->SetEdgeColor(0.23, 0.29, 0.24);
    impl_->actor->GetProperty()->SetLineWidth(1.0);

    impl_->renderer->AddActor(impl_->actor);

    // 数据加载完成后同步工具栏状态，再重建显示管线和默认视角。
    updateScalarChoices();
    toolbar_->setSliceEnabled(sliceEnabled_);
    updateSlicePanelControls();
    toolbar_->setDisplayMode(displayMode_);
    toolbar_->setPerspectiveEnabled(perspectiveEnabled_);
    updateVisualizationPipeline();
    if (resetView) {
        resetCamera();
    } else {
        setPerspectiveEnabled(perspectiveEnabled_);
        vtkCamera* camera = impl_->renderer->GetActiveCamera();
        if (camera && savedCamera) {
            camera->SetPosition(savedPosition);
            camera->SetFocalPoint(savedFocalPoint);
            camera->SetViewUp(savedViewUp);
            camera->SetClippingRange(savedClippingRange);
            camera->SetParallelScale(savedParallelScale);
            camera->SetViewAngle(savedViewAngle);
            camera->SetParallelProjection(savedParallelProjection);
            camera->OrthogonalizeViewUp();
            impl_->window->Render();
        }
    }

    loadedDataPath_ = info.absoluteFilePath();
    loadedDataModified_ = info.lastModified().toMSecsSinceEpoch();
    loadedDataSize_ = info.size();
    topologyPvdPath_ = fromTimeSeries ? loadedPvdPath_ : QString();
    if (!fromTimeSeries) {
        QString message = QStringLiteral("已加载 %1").arg(info.fileName());
        if (!sfmSetLabels.isEmpty()) {
            message += QStringLiteral("；sets：%1")
                           .arg(sfmSetLabels.join(QStringLiteral(", ")));
        }
        emit statusMessage(message);
    }
    return true;
#else
    if (placeholder_) {
        placeholder_->setText(
            QStringLiteral("已选择\n%1\n\n当前构建未找到 VTK，无法渲染。")
                .arg(QFileInfo(path).fileName()));
    }
    if (error) *error = QStringLiteral("当前 sonicGui 构建未启用 VTK");
    Q_UNUSED(resetView)
    Q_UNUSED(fromTimeSeries)
    return false;
#endif
}

bool VTKView::loadPvd(const QString& path, QString* error)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }

    QXmlStreamReader xml(&file);
    QVector<TimeFrame> frames;
    while (!xml.atEnd()) {
        xml.readNext();
        if (!xml.isStartElement()
            || xml.name().compare(QStringLiteral("DataSet"),
                                  Qt::CaseInsensitive) != 0) {
            continue;
        }
        const auto attributes = xml.attributes();
        const QString candidate =
            attributes.value(QStringLiteral("file")).toString();
        if (candidate.isEmpty()) continue;
        const QString dataPath = QFileInfo(candidate).isAbsolute()
            ? candidate : QDir(QFileInfo(path).absolutePath()).filePath(candidate);
        frames.push_back({
            QFileInfo(dataPath).absoluteFilePath(),
            attributes.value(QStringLiteral("timestep")).toString()
        });
    }
    if (xml.hasError()) {
        if (error) {
            *error = QStringLiteral("PVD 解析失败：%1").arg(xml.errorString());
        }
        return false;
    }
    if (frames.isEmpty()) {
        if (error) *error = QStringLiteral("PVD 中没有可显示的时间步");
        return false;
    }

    const bool samePvd = loadedPvdPath_ == QFileInfo(path).absoluteFilePath();
    const int previousCount = timeFrames_.size();
    const QString previousFramePath =
        samePvd && currentFrameIndex_ >= 0 && currentFrameIndex_ < previousCount
            ? timeFrames_.at(currentFrameIndex_).filePath
            : QString();
    const bool wasAtLatest = !samePvd || currentFrameIndex_ < 0
        || currentFrameIndex_ >= previousCount - 1;

    loadedPvdPath_ = QFileInfo(path).absoluteFilePath();
    timeFrames_ = frames;

    int targetIndex = -1;
    if (!wasAtLatest && !previousFramePath.isEmpty()) {
        for (int i = 0; i < timeFrames_.size(); ++i) {
            if (timeFrames_.at(i).filePath == previousFramePath) {
                targetIndex = i;
                break;
            }
        }
    }
    if (targetIndex < 0) targetIndex = timeFrames_.size() - 1;

    const bool loaded =
        loadTimeFrame(targetIndex, error, !samePvd || loadedDataPath_.isEmpty());
    if (loaded && visualizationLoadedHandler_) {
        visualizationLoadedHandler_(loadedPvdPath_);
    }
    return loaded;
}

bool VTKView::loadTimeFrame(int frameIndex, QString* error, bool resetView)
{
    if (frameIndex < 0 || frameIndex >= timeFrames_.size()) {
        if (error) *error = QStringLiteral("时间步索引超出范围");
        return false;
    }

    const TimeFrame& frame = timeFrames_.at(frameIndex);
    const QFileInfo dataInfo(frame.filePath);
    if (!dataInfo.isFile()) {
        if (error) {
            *error = QStringLiteral("PVD 引用的结果尚未生成：%1")
                .arg(frame.filePath);
        }
        return false;
    }

    const bool unchanged = loadedDataPath_ == dataInfo.absoluteFilePath()
        && loadedDataModified_ == dataInfo.lastModified().toMSecsSinceEpoch()
        && loadedDataSize_ == dataInfo.size();
    const bool sameFrame = unchanged && currentFrameIndex_ == frameIndex;
    if (!unchanged) {
        if (!loadDataFile(dataInfo.absoluteFilePath(), error,
                          resetView, true)) {
            return false;
        }
    }

    currentFrameIndex_ = frameIndex;
    updateTimeControls();
    if (!sameFrame) {
        emit statusMessage(
            QStringLiteral("结果时间 %1 · %2")
                .arg(frame.physicalTime.isEmpty()
                         ? QStringLiteral("-") : frame.physicalTime,
                     dataInfo.fileName()));
    }
    return true;
}

void VTKView::watchResultDirectory(const QString& directory)
{
    const QString absolute = directory.isEmpty()
        ? QString() : QFileInfo(directory).absoluteFilePath();
    if (absolute == resultDirectory_) {
        // 工作目录未变时仍安排一次刷新，用于手动切回网格页时拉取最新结果。
        scheduleResultRefresh();
        return;
    }

    // QFileSystemWatcher 不能自动替换路径，切换 case 时先解除旧监听。
    const QStringList watchedFiles = resultWatcher_->files();
    if (!watchedFiles.isEmpty()) resultWatcher_->removePaths(watchedFiles);
    const QStringList watchedDirectories = resultWatcher_->directories();
    if (!watchedDirectories.isEmpty()) {
        resultWatcher_->removePaths(watchedDirectories);
    }

    resultDirectory_ = absolute;
    preferredPvd_.clear();
    watchedPvd_.clear();
    loadedPvdPath_.clear();
    loadedDataPath_.clear();
    loadedDataModified_ = -1;
    loadedDataSize_ = -1;
    topologyPvdPath_.clear();
    timeFrames_.clear();
    currentFrameIndex_ = -1;
    setPlaybackEnabled(false);
    updateTimeControls();
    if (!resultDirectory_.isEmpty() && QDir(resultDirectory_).exists()) {
        resultWatcher_->addPath(resultDirectory_);
    }
    scheduleResultRefresh();
}

void VTKView::setVisualizationLoadedHandler(
    std::function<void(const QString&)> handler)
{
    visualizationLoadedHandler_ = std::move(handler);
}

void VTKView::refreshResult()
{
    if (resultDirectory_.isEmpty()) return;

    QDir directory(resultDirectory_);
    QFileInfo selected;
    const QFileInfo preferred(preferredPvd_);
    if (preferred.isFile()
        && preferred.absolutePath() == directory.absolutePath()) {
        selected = preferred;
    } else {
        const QFileInfoList pvdFiles = directory.entryInfoList(
            {QStringLiteral("*.pvd"), QStringLiteral("*.PVD")},
            QDir::Files,
            QDir::Time);
        for (const QFileInfo& candidate : pvdFiles) {
            // macOS 资源叉文件不是求解结果，必须跳过。
            if (candidate.fileName().startsWith(QStringLiteral("._"))) continue;
            selected = candidate;
            break;
        }
    }
    if (!selected.isFile()) return;

    const QString pvdPath = selected.absoluteFilePath();
    if (watchedPvd_ != pvdPath
        || !resultWatcher_->files().contains(pvdPath)) {
        if (!watchedPvd_.isEmpty()
            && resultWatcher_->files().contains(watchedPvd_)) {
            resultWatcher_->removePath(watchedPvd_);
        }
        watchedPvd_ = pvdPath;
        if (QFileInfo::exists(watchedPvd_)) {
            resultWatcher_->addPath(watchedPvd_);
        }
    }

    QString error;
    if (!loadPvd(pvdPath, &error) && !error.isEmpty()) {
        emit statusMessage(error);
    }
}

void VTKView::activateMeshDefaults()
{
    setPerspectiveEnabled(true);
    if (toolbar_) toolbar_->setViewDirection(ViewerCameraDirection::PositiveZ);
    setCameraDirection(ViewerCameraDirection::PositiveZ);
}

void VTKView::showGeometryScene()
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_ || !impl_->geometryDataReady) return;
    if (impl_->actor) impl_->renderer->RemoveActor(impl_->actor);
    if (impl_->renderer->HasViewProp(impl_->scalarBar.GetPointer())) {
        impl_->renderer->RemoveViewProp(impl_->scalarBar.GetPointer());
    }
    if (!impl_->renderer->HasViewProp(impl_->geometryDefaultFaceActor)) {
        impl_->renderer->AddActor(impl_->geometryDefaultFaceActor);
    }
    if (!impl_->renderer->HasViewProp(impl_->geometryPatchActor)) {
        impl_->renderer->AddActor(impl_->geometryPatchActor);
    }
    if (!impl_->renderer->HasViewProp(impl_->geometryLineActor)) {
        impl_->renderer->AddActor(impl_->geometryLineActor);
    }
    impl_->geometrySceneVisible = true;
    impl_->renderer->ResetCameraClippingRange();
    impl_->window->Render();
#endif
}

void VTKView::showResultScene()
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_) return;
    if (impl_->highlightActor) impl_->renderer->RemoveActor(impl_->highlightActor);
    impl_->highlightActor = nullptr;
    impl_->highlightMapper = nullptr;
    if (impl_->geometryLineActor)
        impl_->renderer->RemoveActor(impl_->geometryLineActor);
    if (impl_->geometryPatchActor)
        impl_->renderer->RemoveActor(impl_->geometryPatchActor);
    if (impl_->geometryDefaultFaceActor)
        impl_->renderer->RemoveActor(impl_->geometryDefaultFaceActor);
    impl_->geometrySceneVisible = false;
    if (impl_->actor && !impl_->renderer->HasViewProp(impl_->actor)) {
        impl_->renderer->AddActor(impl_->actor);
    }
    if (impl_->actor && impl_->surfaceData) updateVisualizationPipeline();
    impl_->renderer->ResetCameraClippingRange();
    impl_->window->Render();
#endif
}

void VTKView::scheduleResultRefresh()
{
    resultReloadTimer_->start();
}

void VTKView::resetCamera()
{
#ifdef SONIC_GUI_WITH_VTK
    setCameraDirection(ViewerCameraDirection::PositiveZ);
#endif
}

void VTKView::setSliceEnabled(bool enabled)
{
    if (enabled) {
        // 打开面板时先同步当前数据范围，避免用户输入超出 bounds。
        updateSlicePanelControls();
        if (slicePanel_) {
            slicePanel_->show();
            slicePanel_->raise();
        }
        return;
    }
    closeSlicePanel();
}

void VTKView::applySliceFromPanel()
{
#ifdef SONIC_GUI_WITH_VTK
    // 坐标框定义切片面经过的点。
    if (sliceXInput_) sliceOrigin_[0] = sliceXInput_->value();
    if (sliceYInput_) sliceOrigin_[1] = sliceYInput_->value();
    if (sliceZInput_) sliceOrigin_[2] = sliceZInput_->value();

    if (impl_ && impl_->renderer) {
        vtkCamera* camera = impl_->renderer->GetActiveCamera();
        if (camera) {
            // 切片法向取当前视线方向：相机焦点 - 相机位置。
            // 用户先切到某个方向再 Slice，就能得到类似 ParaView 的视角切片。
            double position[3];
            double focalPoint[3];
            camera->GetPosition(position);
            camera->GetFocalPoint(focalPoint);
            double normal[3] = {
                focalPoint[0] - position[0],
                focalPoint[1] - position[1],
                focalPoint[2] - position[2]
            };
            if (vtkMath::Normalize(normal) > 1.0e-12) {
                sliceNormal_[0] = normal[0];
                sliceNormal_[1] = normal[1];
                sliceNormal_[2] = normal[2];
            }
        }
    }
#endif
    sliceEnabled_ = true;
    if (slicePanel_) slicePanel_->hide();
    if (toolbar_) toolbar_->setSliceEnabled(true);
    updateSliceData();
    updateVisualizationPipeline();
}

void VTKView::closeSlicePanel()
{
    // 关闭 Slice 后不保留切片视图，直接恢复完整 surface 显示。
    sliceEnabled_ = false;
    if (slicePanel_) slicePanel_->hide();
    if (toolbar_) toolbar_->setSliceEnabled(false);
    updateVisualizationPipeline();
}

void VTKView::setScalarField(const QString& scalarName)
{
    // 只记录目标标量名；真正的 VTK mapper 配置在统一显示管线里完成。
    selectedScalar_ = scalarName;
    updateVisualizationPipeline();
}

void VTKView::setDisplayMode(ViewerDisplayMode mode)
{
    // PhysicalSurface 与 Mesh 共用同一个 surface 数据，
    // 区别只是 actor 是否打开 EdgeVisibility。
    displayMode_ = mode;
    if (toolbar_) toolbar_->setDisplayMode(mode);
    updateVisualizationPipeline();
}

void VTKView::setPerspectiveEnabled(bool enabled)
{
    perspectiveEnabled_ = enabled;
    if (toolbar_) toolbar_->setPerspectiveEnabled(enabled);
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_->renderer) return;
    vtkCamera* camera = impl_->renderer->GetActiveCamera();
    if (!camera) return;
    camera->SetParallelProjection(!enabled);
    if (!enabled && impl_->actor) {
        // 正交投影需要手动设置 ParallelScale，否则模型可能过大或过小。
        double center[3];
        double radius = 1.0;
        boundsCenterAndRadius(impl_->actor, center, &radius);
        camera->SetParallelScale(radius);
    }
    impl_->renderer->ResetCameraClippingRange();
    impl_->window->Render();
#endif
}

void VTKView::setCameraDirection(ViewerCameraDirection direction)
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_->actor) return;
    double center[3];
    double radius = 1.0;
    boundsCenterAndRadius(impl_->actor, center, &radius);

    // axis 是相机从模型中心看过去的方向；up 控制屏幕上方朝向。
    double axis[3] = {1.0, 0.0, 0.0};
    double up[3] = {0.0, 0.0, 1.0};
    switch (direction) {
        case ViewerCameraDirection::PositiveX:
            axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0; break;
        case ViewerCameraDirection::NegativeX:
            axis[0] = -1.0; axis[1] = 0.0; axis[2] = 0.0; break;
        case ViewerCameraDirection::PositiveY:
            axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0; break;
        case ViewerCameraDirection::NegativeY:
            axis[0] = 0.0; axis[1] = -1.0; axis[2] = 0.0; break;
        case ViewerCameraDirection::PositiveZ:
            axis[0] = 0.0; axis[1] = 0.0; axis[2] = 1.0;
            up[0] = 0.0; up[1] = 1.0; up[2] = 0.0;
            break;
        case ViewerCameraDirection::NegativeZ:
            axis[0] = 0.0; axis[1] = 0.0; axis[2] = -1.0;
            up[0] = 0.0; up[1] = 1.0; up[2] = 0.0;
            break;
    }

    vtkCamera* camera = impl_->renderer->GetActiveCamera();
    // 距离取包围球半径的三倍，配合裁剪范围重置可完整显示模型。
    const double distance = radius * 3.0;
    camera->SetFocalPoint(center);
    camera->SetPosition(center[0] + axis[0] * distance,
                        center[1] + axis[1] * distance,
                        center[2] + axis[2] * distance);
    camera->SetViewUp(up);
    camera->OrthogonalizeViewUp();
    camera->SetParallelProjection(!perspectiveEnabled_);
    if (!perspectiveEnabled_) camera->SetParallelScale(radius);
    impl_->renderer->ResetCameraClippingRange();
    impl_->window->Render();
#else
    Q_UNUSED(direction)
#endif
}

void VTKView::setPlaybackEnabled(bool enabled)
{
    if (enabled && timeFrames_.size() > 1) {
        playbackTimer_->start();
        if (toolbar_) toolbar_->setPlaybackEnabled(true);
        return;
    }
    playbackTimer_->stop();
    if (toolbar_) toolbar_->setPlaybackEnabled(false);
}

void VTKView::setPlaybackRate(int framesPerSecond)
{
    playbackRate_ = std::clamp(framesPerSecond, 1, 999);
    playbackTimer_->setInterval(
        std::max(1, qRound(1000.0 / static_cast<double>(playbackRate_))));
}

void VTKView::advancePlayback()
{
    if (timeFrames_.size() <= 1) {
        setPlaybackEnabled(false);
        return;
    }
    const int nextIndex = currentFrameIndex_ >= timeFrames_.size() - 1
        ? 0 : currentFrameIndex_ + 1;
    QString error;
    if (!loadTimeFrame(nextIndex, &error, false)) {
        setPlaybackEnabled(false);
        if (!error.isEmpty()) emit statusMessage(error);
    }
}

void VTKView::showPreviousTimeStep()
{
    if (timeFrames_.isEmpty()) return;
    const int index = currentFrameIndex_ <= 0
        ? 0 : currentFrameIndex_ - 1;
    QString error;
    if (!loadTimeFrame(index, &error, false) && !error.isEmpty()) {
        emit statusMessage(error);
    }
}

void VTKView::showNextTimeStep()
{
    if (timeFrames_.isEmpty()) return;
    const int index = currentFrameIndex_ >= timeFrames_.size() - 1
        ? timeFrames_.size() - 1 : currentFrameIndex_ + 1;
    QString error;
    if (!loadTimeFrame(index, &error, false) && !error.isEmpty()) {
        emit statusMessage(error);
    }
}

void VTKView::showFirstTimeStep()
{
    if (timeFrames_.isEmpty()) return;
    QString error;
    if (!loadTimeFrame(0, &error, false) && !error.isEmpty()) {
        emit statusMessage(error);
    }
}

void VTKView::showLastTimeStep()
{
    if (timeFrames_.isEmpty()) return;
    QString error;
    if (!loadTimeFrame(timeFrames_.size() - 1, &error, false)
        && !error.isEmpty()) {
        emit statusMessage(error);
    }
}

void VTKView::updateTimeControls()
{
    if (!toolbar_) return;
    const QString physicalTime =
        currentFrameIndex_ >= 0 && currentFrameIndex_ < timeFrames_.size()
            ? timeFrames_.at(currentFrameIndex_).physicalTime
            : QString();
    toolbar_->setTimeFrameState(currentFrameIndex_,
                                timeFrames_.size(),
                                physicalTime);
    if (timeFrames_.size() <= 1) setPlaybackEnabled(false);
}

void VTKView::updateVisualizationPipeline()
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_->mapper || !impl_->actor || !impl_->surfaceData) return;
    // 显示数据在“完整 surface”和“当前切片”之间切换。
    vtkPolyData* data = sliceEnabled_ && impl_->sliceData
        ? impl_->sliceData.GetPointer()
        : impl_->surfaceData.GetPointer();
    impl_->displayData = data;
    impl_->mapper->SetInputData(data);

    vtkDataArray* scalar = findScalarArray(data, selectedScalar_);
    if (!scalar) {
        // STL 或没有目标标量的数据只能用固定颜色显示。
        impl_->mapper->ScalarVisibilityOff();
        impl_->actor->GetProperty()->SetColor(0.42, 0.53, 0.44);
        if (impl_->renderer->HasViewProp(impl_->scalarBar.GetPointer())) {
            impl_->renderer->RemoveViewProp(impl_->scalarBar.GetPointer());
        }
    } else {
        // 标量着色使用点数据数组，范围随当前显示数据（整面或切片）更新。
        double range[2];
        scalar->GetRange(range);
        impl_->lookupTable->SetRange(range);
        impl_->lookupTable->Build();
        impl_->mapper->SetLookupTable(impl_->lookupTable);
        impl_->mapper->SetScalarModeToUsePointFieldData();
        impl_->mapper->SelectColorArray(scalar->GetName());
        impl_->mapper->SetScalarRange(range);
        impl_->mapper->ScalarVisibilityOn();
        impl_->scalarBar->SetTitle(scalar->GetName());
        impl_->scalarBar->SetLookupTable(impl_->lookupTable);
        if (!impl_->renderer->HasViewProp(impl_->scalarBar.GetPointer())) {
            impl_->renderer->AddViewProp(impl_->scalarBar.GetPointer());
        }
    }

    const bool meshMode = displayMode_ == ViewerDisplayMode::Mesh;
    impl_->actor->GetProperty()->SetRepresentationToSurface();
    if (meshMode) {
        // ParaView 的 Surface With Edges：仍然渲染物理面，只额外显示单元边。
        impl_->actor->GetProperty()->EdgeVisibilityOn();
        impl_->actor->GetProperty()->SetEdgeColor(0.07, 0.09, 0.08);
        impl_->actor->GetProperty()->SetLineWidth(0.8);
    } else {
        impl_->actor->GetProperty()->EdgeVisibilityOff();
        impl_->actor->GetProperty()->SetLineWidth(1.0);
    }
    impl_->window->Render();
#endif
}

void VTKView::updateScalarChoices()
{
#ifdef SONIC_GUI_WITH_VTK
    // 下拉框选项来自完整 surface。切片数据只是显示子集，不单独改变列表。
    const QStringList names = scalarNames(impl_->surfaceData);
    selectedScalar_ = chooseScalar(names, selectedScalar_);
    toolbar_->setScalarChoices(names, selectedScalar_);
#endif
}

void VTKView::updateSlicePanelControls()
{
#ifdef SONIC_GUI_WITH_VTK
    // 更新输入框时屏蔽 valueChanged，避免刷新范围时误触发外部逻辑。
    auto syncInput = [](QDoubleSpinBox* input,
                        double minimum,
                        double maximum,
                        double value) {
        if (!input) return;
        const QSignalBlocker blocker(input);
        input->setRange(minimum, maximum);
        input->setSingleStep(spinStep(minimum, maximum));
        input->setValue(value);
    };
    syncInput(sliceXInput_, sliceBounds_[0], sliceBounds_[1], sliceOrigin_[0]);
    syncInput(sliceYInput_, sliceBounds_[2], sliceBounds_[3], sliceOrigin_[1]);
    syncInput(sliceZInput_, sliceBounds_[4], sliceBounds_[5], sliceOrigin_[2]);
#endif
}

void VTKView::updateSliceData()
{
#ifdef SONIC_GUI_WITH_VTK
    // 切片缓存独立保存；是否显示它由 updateVisualizationPipeline 决定。
    impl_->sliceData = sliceFromDataSet(
        impl_->sourceData, sliceOrigin_, sliceNormal_);
#endif
}


void VTKView::highlightBlock(int blockId)
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_ || !impl_->geometryBlockData) return;
    clearHighlight();
    const auto extracted = cellsWithId(impl_->geometryBlockData, "blockID", blockId);
    if (extracted->GetNumberOfCells() == 0) return;
    impl_->highlightMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    impl_->highlightMapper->SetInputData(extracted);
    impl_->highlightMapper->ScalarVisibilityOff();
    impl_->highlightActor = vtkSmartPointer<vtkActor>::New();
    impl_->highlightActor->SetMapper(impl_->highlightMapper);
    impl_->highlightActor->GetProperty()->SetColor(1.0, 0.0, 0.0);      // 红色
    impl_->highlightActor->GetProperty()->SetLineWidth(2.0);
    impl_->highlightActor->GetProperty()->SetOpacity(1.0);
    impl_->highlightActor->GetProperty()->SetAmbient(1.0);
    impl_->highlightActor->GetProperty()->SetDiffuse(0.0);
    impl_->highlightActor->GetProperty()->EdgeVisibilityOn();
    impl_->highlightActor->GetProperty()->SetEdgeColor(0.05, 0.05, 0.05);
    impl_->renderer->AddActor(impl_->highlightActor);
    if (impl_->geometryLineActor) { impl_->geometryLineActor->GetProperty()->SetColor(0.05,0.05,0.05); impl_->geometryLineActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryPatchActor) { impl_->geometryPatchActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryPatchActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryDefaultFaceActor) { impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(0.8); }
    impl_->highlightActive = true;
    impl_->widget->renderWindow()->Render();
#endif
}

void VTKView::highlightPatch(int patchId)
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_ || !impl_->geometryPatchData) return;
    clearHighlight();
    const auto extracted = cellsWithId(impl_->geometryPatchData, "patchID", patchId);
    if (extracted->GetNumberOfCells() == 0) return;
    impl_->highlightMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    impl_->highlightMapper->SetInputData(extracted);
    impl_->highlightMapper->ScalarVisibilityOff();
    impl_->highlightActor = vtkSmartPointer<vtkActor>::New();
    impl_->highlightActor->SetMapper(impl_->highlightMapper);
    impl_->highlightActor->GetProperty()->SetColor(0.0, 0.6, 0.0);       // 正面绿色
    impl_->highlightActor->GetProperty()->SetOpacity(0.9);
    impl_->highlightActor->GetProperty()->SetAmbient(1.0);
    impl_->highlightActor->GetProperty()->SetDiffuse(0.0);
    impl_->highlightActor->GetProperty()->SetEdgeVisibility(true);
    impl_->highlightActor->GetProperty()->SetEdgeColor(0.2, 0.2, 0.2);
    auto backface = vtkSmartPointer<vtkProperty>::New();
    backface->SetColor(0.0, 0.0, 0.8); // 反面蓝色
    backface->SetOpacity(0.9);
    backface->SetAmbient(1.0);
    backface->SetDiffuse(0.0);
    impl_->highlightActor->SetBackfaceProperty(backface);
    impl_->renderer->AddActor(impl_->highlightActor);
    if (impl_->geometryLineActor) { impl_->geometryLineActor->GetProperty()->SetColor(0.05,0.05,0.05); impl_->geometryLineActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryPatchActor) { impl_->geometryPatchActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryPatchActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryDefaultFaceActor) { impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(0.8); }
    impl_->highlightActive = true;
    impl_->widget->renderWindow()->Render();
#endif
}

void VTKView::highlightFace(int patchId, int faceId)
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_ || !impl_->geometryPatchData) return;
    clearHighlight();
    const auto patchOnly = cellsWithId(impl_->geometryPatchData, "patchID", patchId);
    if (patchOnly->GetNumberOfCells() == 0) return;
    const auto extracted = cellsWithId(patchOnly, "faceID", faceId);
    if (extracted->GetNumberOfCells() == 0) return;
    impl_->highlightMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    impl_->highlightMapper->SetInputData(extracted);
    impl_->highlightMapper->ScalarVisibilityOff();
    impl_->highlightActor = vtkSmartPointer<vtkActor>::New();
    impl_->highlightActor->SetMapper(impl_->highlightMapper);
    impl_->highlightActor->GetProperty()->SetColor(0.6, 0.0, 0.8);       // 正面紫色
    impl_->highlightActor->GetProperty()->SetOpacity(0.95);
    impl_->highlightActor->GetProperty()->SetAmbient(1.0);
    impl_->highlightActor->GetProperty()->SetDiffuse(0.0);
    impl_->highlightActor->GetProperty()->SetEdgeVisibility(true);
    impl_->highlightActor->GetProperty()->SetEdgeColor(0.1, 0.1, 0.1);
    auto backface = vtkSmartPointer<vtkProperty>::New();
    backface->SetColor(0.6, 0.3, 0.0); // 反面棕色
    backface->SetOpacity(0.95);
    backface->SetAmbient(1.0);
    backface->SetDiffuse(0.0);
    impl_->highlightActor->SetBackfaceProperty(backface);
    impl_->renderer->AddActor(impl_->highlightActor);
    if (impl_->geometryLineActor) { impl_->geometryLineActor->GetProperty()->SetColor(0.05,0.05,0.05); impl_->geometryLineActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryPatchActor) { impl_->geometryPatchActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryPatchActor->GetProperty()->SetOpacity(0.8); }
    if (impl_->geometryDefaultFaceActor) { impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0,1.0,1.0); impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(0.8); }
    impl_->highlightActive = true;
    impl_->widget->renderWindow()->Render();
#endif
}

void VTKView::clearHighlight()
{
#ifdef SONIC_GUI_WITH_VTK
    if (!impl_) return;
    if (impl_->highlightActor) {
        impl_->renderer->RemoveActor(impl_->highlightActor);
        impl_->highlightActor = nullptr;
    }
    impl_->highlightMapper = nullptr;
    if (impl_->geometryLineActor) { impl_->geometryLineActor->GetProperty()->SetColor(0.05,0.05,0.05); impl_->geometryLineActor->GetProperty()->SetOpacity(1.0); }
    if (impl_->geometryPatchActor) {
        impl_->geometryPatchActor->GetProperty()->SetColor(1.0,1.0,1.0);
        impl_->geometryPatchActor->GetProperty()->SetOpacity(1.0);
        auto* back = impl_->geometryPatchActor->GetBackfaceProperty();
        if (back) { back->SetColor(1.0,1.0,1.0); back->SetOpacity(1.0); }
    }
    if (impl_->geometryDefaultFaceActor) {
        impl_->geometryDefaultFaceActor->GetProperty()->SetColor(1.0,1.0,1.0);
        impl_->geometryDefaultFaceActor->GetProperty()->SetOpacity(1.0);
        auto* back = impl_->geometryDefaultFaceActor->GetBackfaceProperty();
        if (back) { back->SetColor(1.0,1.0,1.0); back->SetOpacity(1.0); }
    }
    impl_->highlightActive = false;
    impl_->widget->renderWindow()->Render();
#endif
}
}
