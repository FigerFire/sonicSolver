/// @file SF_editor.cpp
/// @brief GUI 几何编辑与参数面板实现。

#include "SF_editor.h"

#include <QFile>
#include <QFileInfo>
#include <QGraphicsPathItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QHash>
#include <QHBoxLayout>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPushButton>
#include <QMap>
#include <QRegularExpression>
#include <QResizeEvent>
#include <QTimer>
#include <QToolButton>
#include <QVBoxLayout>

#include <cmath>

namespace SF::GUI {
namespace {

struct MeshVertex {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct DisplayEdge {
    int start = -1;
    int end = -1;
};

struct DisplayArc {
    int start = -1;
    int end = -1;
    QPointF middle;
};

struct DisplayBlock {
    QVector<int> vertices;
    QColor color;
    QString label;
};

struct DisplayPatch {
    QString name;
    QColor color;
    QVector<DisplayEdge> edges;
};

struct BlockMeshDisplay {
    QVector<MeshVertex> vertices;
    QVector<DisplayEdge> edges;
    QVector<DisplayArc> arcs;
    QVector<DisplayBlock> blocks;
    QVector<DisplayPatch> patches;
};

bool extractArrayValue(const QString& text,
                       const QString& key,
                       QString* value)
{
    const QRegularExpression assignment(
        QStringLiteral(R"(\b%1\s*=)")
            .arg(QRegularExpression::escape(key)));
    const QRegularExpressionMatch match = assignment.match(text);
    if (!match.hasMatch()) return false;

    const int open = text.indexOf(QLatin1Char('['), match.capturedEnd());
    if (open < 0) return false;

    int depth = 0;
    bool quoted = false;
    QChar quote;
    for (int i = open; i < text.size(); ++i) {
        const QChar ch = text.at(i);
        if (ch == QLatin1Char('"') || ch == QLatin1Char('\'')) {
            if (!quoted) {
                quoted = true;
                quote = ch;
            } else if (quote == ch
                       && (i == 0 || text.at(i - 1) != QLatin1Char('\\'))) {
                quoted = false;
            }
            continue;
        }
        if (quoted) continue;
        if (ch == QLatin1Char('[')) ++depth;
        if (ch != QLatin1Char(']')) continue;
        --depth;
        if (depth == 0) {
            *value = text.mid(open, i - open + 1);
            return true;
        }
    }
    return false;
}

QRegularExpression numberPattern()
{
    return QRegularExpression(
        QStringLiteral(R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"));
}

QString stripFoamComments(QString text)
{
    text.remove(QRegularExpression(QStringLiteral(R"(/\*[\s\S]*?\*/)")));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\n]*)")));
    return text;
}

bool isWordChar(QChar ch)
{
    return ch.isLetterOrNumber() || ch == QLatin1Char('_')
        || ch == QLatin1Char('-') || ch == QLatin1Char('.');
}

int findMatching(const QString& text, int openPos, QChar open, QChar close)
{
    int depth = 0;
    bool quoted = false;
    QChar quote;
    for (int i = openPos; i < text.size(); ++i) {
        const QChar ch = text.at(i);
        if (ch == QLatin1Char('"') || ch == QLatin1Char('\'')) {
            if (!quoted) {
                quoted = true;
                quote = ch;
            } else if (quote == ch
                       && (i == 0 || text.at(i - 1) != QLatin1Char('\\'))) {
                quoted = false;
            }
            continue;
        }
        if (quoted) continue;
        if (ch == open) ++depth;
        else if (ch == close) {
            --depth;
            if (depth == 0) return i;
        }
    }
    return -1;
}

QString foamBlockAfterKey(const QString& text,
                          const QString& key,
                          QChar open,
                          QChar close)
{
    const QRegularExpression keyPattern(
        QStringLiteral(R"(\b%1\b)").arg(QRegularExpression::escape(key)));
    int searchFrom = 0;
    while (true) {
        const QRegularExpressionMatch match = keyPattern.match(text, searchFrom);
        if (!match.hasMatch()) return {};
        int p = match.capturedEnd();
        while (p < text.size() && text.at(p).isSpace()) ++p;
        if (p < text.size() && text.at(p) == open) {
            const int end = findMatching(text, p, open, close);
            if (end >= 0) return text.mid(p + 1, end - p - 1);
        }
        searchFrom = match.capturedEnd();
    }
}

QVector<int> parseIndexList(const QString& text)
{
    QVector<int> values;
    const QRegularExpression indexPattern(QStringLiteral(R"(\d+)"));
    auto matches = indexPattern.globalMatch(text);
    while (matches.hasNext()) {
        values.push_back(matches.next().captured(0).toInt());
    }
    return values;
}

QPointF projectVertex(const MeshVertex& vertex)
{
    return QPointF(vertex.x, -vertex.y);
}

QPainterPath sampledArcPath(const QPointF& a,
                            const QPointF& b,
                            const QPointF& c)
{
    const double ax = a.x();
    const double ay = a.y();
    const double bx = b.x();
    const double by = b.y();
    const double cx = c.x();
    const double cy = c.y();
    const double d =
        2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

    QPainterPath path(a);
    if (std::abs(d) < 1.0e-12) {
        path.quadTo(b, c);
        return path;
    }

    const double a2 = ax * ax + ay * ay;
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    const QPointF center(
        (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d,
        (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d);
    const double radius = std::hypot(ax - center.x(), ay - center.y());
    if (radius <= 0.0) return path;

    auto angle = [&center](const QPointF& p) {
        return std::atan2(p.y() - center.y(), p.x() - center.x());
    };
    double start = angle(a);
    double mid = angle(b);
    double end = angle(c);
    auto normalizePositive = [](double x) {
        constexpr double twoPi = 6.28318530717958647692;
        while (x < 0.0) x += twoPi;
        while (x >= twoPi) x -= twoPi;
        return x;
    };
    auto ccwDelta = [&](double from, double to) {
        return normalizePositive(to - from);
    };
    const double ccwEnd = ccwDelta(start, end);
    const double ccwMid = ccwDelta(start, mid);
    constexpr double twoPi = 6.28318530717958647692;
    double sweep = ccwMid <= ccwEnd ? ccwEnd : ccwEnd - twoPi;
    const int samples = 32;
    for (int i = 1; i <= samples; ++i) {
        const double t = static_cast<double>(i) / samples;
        const double theta = start + sweep * t;
        path.lineTo(center.x() + radius * std::cos(theta),
                    center.y() + radius * std::sin(theta));
    }
    return path;
}

QVector<DisplayEdge> faceEdges(const QVector<int>& face)
{
    QVector<DisplayEdge> edges;
    if (face.size() < 2) return edges;
    for (int i = 0; i < face.size(); ++i) {
        edges.push_back({face.at(i), face.at((i + 1) % face.size())});
    }
    return edges;
}

QColor indexedColor(int index)
{
    static const QVector<QColor> colors{
        QColor(255, 204, 0),
        QColor(52, 199, 89),
        QColor(0, 122, 255),
        QColor(255, 149, 0),
        QColor(175, 82, 222),
        QColor(255, 45, 85),
        QColor(90, 200, 250),
        QColor(162, 132, 94)
    };
    return colors.at(index % colors.size());
}

bool parseBlockMeshDict(const QString& rawText,
                        BlockMeshDisplay* mesh,
                        QString* error)
{
    const QString text = stripFoamComments(rawText);
    const QString number = numberPattern().pattern();
    const QRegularExpression pointPattern(
        QStringLiteral(R"(\(\s*(%1)\s+(%1)\s+(%1)\s*\))").arg(number));

    const QString verticesBody =
        foamBlockAfterKey(text, QStringLiteral("vertices"),
                          QLatin1Char('('), QLatin1Char(')'));
    if (verticesBody.isEmpty()) {
        if (error) *error = QStringLiteral("blockMeshDict 未找到 vertices");
        return false;
    }

    auto pointMatches = pointPattern.globalMatch(verticesBody);
    while (pointMatches.hasNext()) {
        const QRegularExpressionMatch match = pointMatches.next();
        mesh->vertices.push_back({
            match.captured(1).toDouble(),
            match.captured(2).toDouble(),
            match.captured(3).toDouble()
        });
    }
    if (mesh->vertices.isEmpty()) {
        if (error) *error = QStringLiteral("blockMeshDict vertices 为空");
        return false;
    }

    const QString blocksBody =
        foamBlockAfterKey(text, QStringLiteral("blocks"),
                          QLatin1Char('('), QLatin1Char(')'));
    const QRegularExpression hexPattern(
        QStringLiteral(R"(\bhex\s*\(([^\)]*)\))"),
        QRegularExpression::CaseInsensitiveOption);
    auto blockMatches = hexPattern.globalMatch(blocksBody);
    int blockIndex = 0;
    while (blockMatches.hasNext()) {
        const QVector<int> indices =
            parseIndexList(blockMatches.next().captured(1));
        if (indices.size() < 2) continue;
        DisplayBlock block;
        block.vertices = indices;
        block.color = indexedColor(blockIndex);
        block.label = QStringLiteral("block %1").arg(blockIndex);
        mesh->blocks.push_back(block);
        ++blockIndex;
    }

    const QString edgesBody =
        foamBlockAfterKey(text, QStringLiteral("edges"),
                          QLatin1Char('('), QLatin1Char(')'));
    const QRegularExpression arcPattern(
        QStringLiteral(
            R"(\barc\s+(\d+)\s+(\d+)\s*\(\s*(%1)\s+(%1)\s+(%1)\s*\))")
            .arg(number),
        QRegularExpression::CaseInsensitiveOption);
    auto arcMatches = arcPattern.globalMatch(edgesBody);
    while (arcMatches.hasNext()) {
        const QRegularExpressionMatch match = arcMatches.next();
        const int start = match.captured(1).toInt();
        const int end = match.captured(2).toInt();
        mesh->arcs.push_back({
            start,
            end,
            QPointF(match.captured(3).toDouble(),
                    -match.captured(4).toDouble())
        });
    }

    const QString boundaryBody =
        foamBlockAfterKey(text, QStringLiteral("boundary"),
                          QLatin1Char('('), QLatin1Char(')'));
    int searchFrom = 0;
    int patchIndex = 0;
    const QRegularExpression patchStart(QStringLiteral(R"(\b([A-Za-z_][A-Za-z0-9_.-]*)\s*\{)"));
    while (true) {
        const QRegularExpressionMatch match =
            patchStart.match(boundaryBody, searchFrom);
        if (!match.hasMatch()) break;
        const int bodyOpen = match.capturedEnd(0) - 1;
        const int bodyClose =
            findMatching(boundaryBody, bodyOpen,
                         QLatin1Char('{'), QLatin1Char('}'));
        if (bodyClose < 0) break;

        const QString patchBody =
            boundaryBody.mid(bodyOpen + 1, bodyClose - bodyOpen - 1);
        const QString facesBody =
            foamBlockAfterKey(patchBody, QStringLiteral("faces"),
                              QLatin1Char('('), QLatin1Char(')'));
        if (!facesBody.isEmpty()) {
            DisplayPatch patch;
            patch.name = match.captured(1);
            patch.color = indexedColor(patchIndex);
            const QRegularExpression facePattern(QStringLiteral(R"(\(([0-9\s]+)\))"));
            auto faceMatches = facePattern.globalMatch(facesBody);
            while (faceMatches.hasNext()) {
                const QVector<int> face =
                    parseIndexList(faceMatches.next().captured(1));
                patch.edges += faceEdges(face);
            }
            if (!patch.edges.isEmpty()) {
                mesh->patches.push_back(patch);
                ++patchIndex;
            }
        }
        searchFrom = bodyClose + 1;
    }

    return true;
}

} // namespace

class GeometryCanvas final : public QGraphicsView {
public:
    enum class Mode { Select, Line, Arc };

    explicit GeometryCanvas(QWidget* parent = nullptr)
        : QGraphicsView(parent)
    {
        setScene(&scene_);
        scene_.setSceneRect(-450, -225, 900, 450);
        setRenderHint(QPainter::Antialiasing);
        setDragMode(QGraphicsView::RubberBandDrag);
        setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
        setBackgroundBrush(QColor(239, 232, 220));
        setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    }

    void setMode(Mode mode)
    {
        mode_ = mode;
        pending_.clear();
        setDragMode(mode == Mode::Select
                        ? QGraphicsView::RubberBandDrag
                        : QGraphicsView::NoDrag);
    }

    void clearSketch()
    {
        scene_.clear();
        points_.clear();
        arcs_.clear();
        pending_.clear();
        drawGrid();
    }

    GeometrySketch sketch() const
    {
        return {points_, arcs_};
    }

    void setSketch(const GeometrySketch& sketch)
    {
        scene_.clear();
        pending_.clear();
        points_ = sketch.points;
        arcs_ = sketch.arcs;

        for (int i = 0; i < points_.size(); ++i) {
            addEndpoint(points_.at(i), i);
        }

        auto hasArc = [this](int start, int end) {
            for (const GeometryArc& arc : arcs_) {
                if ((arc.start == start && arc.end == end)
                    || (arc.start == end && arc.end == start)) {
                    return true;
                }
            }
            return false;
        };
        for (int i = 0; i < points_.size(); ++i) {
            const int next = (i + 1) % points_.size();
            if (points_.size() > 1 && !hasArc(i, next)) {
                QPen pen(QColor(74, 164, 255), 2.0);
                pen.setCosmetic(true);
                scene_.addLine(QLineF(points_.at(i), points_.at(next)),
                               pen);
            }
        }
        for (const GeometryArc& arc : arcs_) {
            if (arc.start < 0 || arc.start >= points_.size()
                || arc.end < 0 || arc.end >= points_.size()) {
                continue;
            }
            QPainterPath path(points_.at(arc.start));
            path.quadTo(arc.middle, points_.at(arc.end));
            QPen pen(QColor(255, 159, 10), 2.0);
            pen.setCosmetic(true);
            scene_.addPath(path, pen);
        }
        QTimer::singleShot(0, this, [this] { fitSketch(); });
    }

    void setBlockMesh(const BlockMeshDisplay& mesh)
    {
        scene_.clear();
        pending_.clear();
        points_.clear();
        arcs_.clear();

        for (const MeshVertex& vertex : mesh.vertices) {
            points_.push_back(projectVertex(vertex));
        }

        for (const DisplayBlock& block : mesh.blocks) {
            QPen pen(block.color, 5.0);
            pen.setCosmetic(true);
            pen.setCapStyle(Qt::RoundCap);
            for (const DisplayEdge& edge : blockEdges(block.vertices)) {
                drawEdge(edge, pen);
            }
        }

        QPen basePen(QColor(74, 164, 255), 1.6);
        basePen.setCosmetic(true);
        for (const DisplayEdge& edge : mesh.edges) drawEdge(edge, basePen);

        for (const DisplayArc& arc : mesh.arcs) {
            if (!validIndex(arc.start) || !validIndex(arc.end)) continue;
            QPen pen(QColor(255, 159, 10), 2.4);
            pen.setCosmetic(true);
            scene_.addPath(sampledArcPath(points_.at(arc.start),
                                          arc.middle,
                                          points_.at(arc.end)),
                           pen);
        }

        for (const DisplayPatch& patch : mesh.patches) {
            QPen pen(patch.color, 3.2);
            pen.setCosmetic(true);
            pen.setCapStyle(Qt::RoundCap);
            for (const DisplayEdge& edge : patch.edges) drawEdge(edge, pen);
        }

        for (int i = 0; i < mesh.vertices.size(); ++i) {
            addEndpoint(projectVertex(mesh.vertices.at(i)), i);
        }
        drawLegend(mesh);
        QTimer::singleShot(0, this, [this] { fitSketch(); });
    }

protected:
    void drawBackground(QPainter* painter, const QRectF& rect) override
    {
        QGraphicsView::drawBackground(painter, rect);
        const int grid = 25;
        const int left = static_cast<int>(std::floor(rect.left() / grid)) * grid;
        const int top = static_cast<int>(std::floor(rect.top() / grid)) * grid;
        QPen pen(QColor(205, 195, 178));
        pen.setWidthF(0.0);
        painter->setPen(pen);
        for (int x = left; x < rect.right(); x += grid) {
            painter->drawLine(QLineF(x, rect.top(), x, rect.bottom()));
        }
        for (int y = top; y < rect.bottom(); y += grid) {
            painter->drawLine(QLineF(rect.left(), y, rect.right(), y));
        }
    }

    void mousePressEvent(QMouseEvent* event) override
    {
        if (event->button() != Qt::LeftButton || mode_ == Mode::Select) {
            QGraphicsView::mousePressEvent(event);
            return;
        }

        const QPointF point = mapToScene(event->pos());
        pending_.push_back(point);
        addMarker(point, pending_.size());

        if (mode_ == Mode::Line && pending_.size() == 2) {
            const int start = endpointIndex(pending_.at(0));
            const int end = endpointIndex(pending_.at(1));
            QPen pen(QColor(74, 164, 255), 2.0);
            pen.setCosmetic(true);
            scene_.addLine(QLineF(points_.at(start), points_.at(end)),
                           pen);
            pending_.clear();
        } else if (mode_ == Mode::Arc && pending_.size() == 3) {
            const int start = endpointIndex(pending_.at(0));
            const int end = endpointIndex(pending_.at(2));
            QPainterPath path(points_.at(start));
            path.quadTo(pending_.at(1), points_.at(end));
            QPen pen(QColor(255, 159, 10), 2.0);
            pen.setCosmetic(true);
            scene_.addPath(path, pen);
            arcs_.push_back({start, end, pending_.at(1)});
            pending_.clear();
        }
    }

    void wheelEvent(QWheelEvent* event) override
    {
        const double factor = event->angleDelta().y() > 0 ? 1.12 : 0.89;
        scale(factor, factor);
    }

    void resizeEvent(QResizeEvent* event) override
    {
        QGraphicsView::resizeEvent(event);
        if (!points_.isEmpty()) {
            QTimer::singleShot(0, this, [this] { fitSketch(); });
        }
    }

private:
    void fitSketch()
    {
        if (points_.isEmpty()) return;
        qreal minX = points_.first().x();
        qreal maxX = minX;
        qreal minY = points_.first().y();
        qreal maxY = minY;
        for (const QPointF& point : points_) {
            minX = std::min(minX, point.x());
            maxX = std::max(maxX, point.x());
            minY = std::min(minY, point.y());
            maxY = std::max(maxY, point.y());
        }
        QRectF bounds(QPointF(minX, minY), QPointF(maxX, maxY));
        const qreal margin =
            std::max<qreal>(1.0, std::max(bounds.width(), bounds.height())
                                  * 0.15);
        bounds.adjust(-margin, -margin, margin, margin);
        if (!bounds.isEmpty()) fitInView(bounds, Qt::KeepAspectRatio);
    }

    bool validIndex(int index) const
    {
        return index >= 0 && index < points_.size();
    }

    QVector<DisplayEdge> blockEdges(const QVector<int>& vertices) const
    {
        if (vertices.size() != 8) return faceEdges(vertices);
        static const int pairs[12][2] = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0},
            {4, 5}, {5, 6}, {6, 7}, {7, 4},
            {0, 4}, {1, 5}, {2, 6}, {3, 7}
        };
        QVector<DisplayEdge> edges;
        for (const auto& pair : pairs) {
            edges.push_back({vertices.at(pair[0]), vertices.at(pair[1])});
        }
        return edges;
    }

    void drawEdge(const DisplayEdge& edge, const QPen& pen)
    {
        if (!validIndex(edge.start) || !validIndex(edge.end)) return;
        scene_.addLine(QLineF(points_.at(edge.start), points_.at(edge.end)),
                       pen);
    }

    int endpointIndex(const QPointF& point)
    {
        constexpr double snapDistance = 10.0;
        for (int i = 0; i < points_.size(); ++i) {
            const QPointF delta = points_.at(i) - point;
            if (std::hypot(delta.x(), delta.y()) <= snapDistance) return i;
        }
        points_.push_back(point);
        addEndpoint(point, points_.size() - 1);
        return points_.size() - 1;
    }

    void addEndpoint(const QPointF& point, int index)
    {
        constexpr double radius = 4.5;
        auto* marker = scene_.addEllipse(-radius,
                                         -radius,
                                         radius * 2,
                                         radius * 2,
                                         QPen(Qt::NoPen),
                                         QBrush(QColor(10, 132, 255)));
        marker->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        marker->setPos(point);
        marker->setToolTip(QStringLiteral("顶点 %1").arg(index));

        auto* label = scene_.addText(QString::number(index));
        label->setDefaultTextColor(QColor(38, 38, 38));
        label->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        label->setPos(point + QPointF(5.0, -18.0));
    }

    void addMarker(const QPointF& point, int step)
    {
        constexpr double radius = 2.5;
        auto* marker = scene_.addEllipse(-radius,
                                         -radius,
                                         radius * 2,
                                         radius * 2,
                                         QPen(Qt::NoPen),
                                         QBrush(QColor(255, 214, 10)));
        marker->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        marker->setPos(point);
        marker->setToolTip(QStringLiteral("当前步骤 %1").arg(step));
    }

    void drawGrid() {}

    void drawLegend(const BlockMeshDisplay& mesh)
    {
        if (mesh.patches.isEmpty()) return;
        qreal minX = 0.0;
        qreal minY = 0.0;
        qreal maxY = 0.0;
        if (!points_.isEmpty()) {
            minX = points_.first().x();
            minY = points_.first().y();
            maxY = minY;
            for (const QPointF& point : points_) {
                minX = std::min(minX, point.x());
                minY = std::min(minY, point.y());
                maxY = std::max(maxY, point.y());
            }
        }
        const QPointF origin(minX, maxY + 36.0);
        for (int i = 0; i < mesh.patches.size(); ++i) {
            const DisplayPatch& patch = mesh.patches.at(i);
            const QPointF p = origin + QPointF(0.0, i * 22.0);
            auto* swatch = scene_.addRect(
                QRectF(p.x(), p.y(), 18.0, 8.0),
                QPen(Qt::NoPen), QBrush(patch.color));
            swatch->setFlag(QGraphicsItem::ItemIgnoresTransformations);
            auto* text = scene_.addText(
                QStringLiteral("%1 = %2").arg(i + 1).arg(patch.name));
            text->setDefaultTextColor(QColor(48, 48, 48));
            text->setFlag(QGraphicsItem::ItemIgnoresTransformations);
            text->setPos(p + QPointF(24.0, -8.0));
        }
    }

    QGraphicsScene scene_;
    Mode mode_ = Mode::Select;
    QVector<QPointF> points_;
    QVector<GeometryArc> arcs_;
    QVector<QPointF> pending_;
};

GeometryEditor::GeometryEditor(QWidget* parent)
    : QWidget(parent)
{
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(0, 0, 0, 0);
    root->setSpacing(0);

    auto* toolbar = new QWidget(this);
    toolbar->setObjectName(QStringLiteral("viewerToolbar"));
    auto* row = new QHBoxLayout(toolbar);
    row->setContentsMargins(12, 7, 12, 7);

    auto* title = new QLabel(QStringLiteral("几何草图"), toolbar);
    title->setObjectName(QStringLiteral("panelTitle"));
    row->addWidget(title);
    row->addSpacing(12);

    row->addStretch();

    auto* clearButton = new QToolButton(toolbar);
    clearButton->setText(QStringLiteral("清空"));
    auto* saveButton = new QPushButton(QStringLiteral("生成单块参数"), toolbar);
    saveButton->setObjectName(QStringLiteral("primaryButton"));
    row->addWidget(clearButton);
    row->addWidget(saveButton);
    root->addWidget(toolbar);

    canvas_ = new GeometryCanvas(this);
    root->addWidget(canvas_, 1);

    statusLabel_ = new QLabel(
        QStringLiteral("绘制四个边界顶点；圆弧按起点、过点、终点依次点击。"),
        this);
    statusLabel_->setObjectName(QStringLiteral("mutedLabel"));
    statusLabel_->setContentsMargins(12, 7, 12, 9);
    statusLabel_->setMinimumHeight(32);
    root->addWidget(statusLabel_);

    connect(clearButton, &QToolButton::clicked,
            canvas_, &GeometryCanvas::clearSketch);
    connect(saveButton, &QPushButton::clicked, this, [this] {
        const GeometrySketch value = canvas_->sketch();
        statusLabel_->setText(
            QStringLiteral("当前草图：%1 个边界顶点，%2 条圆弧")
                .arg(value.points.size())
                .arg(value.arcs.size()));
        emit saveRequested(value);
    });
}

void GeometryEditor::setSelectMode()
{
    canvas_->setMode(GeometryCanvas::Mode::Select);
}

void GeometryEditor::setLineMode()
{
    canvas_->setMode(GeometryCanvas::Mode::Line);
}

void GeometryEditor::setArcMode()
{
    canvas_->setMode(GeometryCanvas::Mode::Arc);
}

bool GeometryEditor::loadFile(const QString& path, QString* error)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }
    const QString text = QString::fromUtf8(file.readAll());

    const QFileInfo info(path);
    if (info.fileName().compare(QStringLiteral("blockMeshDict"),
                                Qt::CaseInsensitive) == 0
        || text.contains(QStringLiteral("object      blockMeshDict"),
                         Qt::CaseInsensitive)
        || text.contains(QStringLiteral("object blockMeshDict"),
                         Qt::CaseInsensitive)) {
        BlockMeshDisplay mesh;
        if (!parseBlockMeshDict(text, &mesh, error)) return false;
        canvas_->setBlockMesh(mesh);
        statusLabel_->setText(
            QStringLiteral("已读取 blockMeshDict：%1 个顶点，%2 个 block，%3 个 patch")
                .arg(mesh.vertices.size())
                .arg(mesh.blocks.size())
                .arg(mesh.patches.size()));
        return true;
    }

    const QString number = numberPattern().pattern();
    const QRegularExpression pointPattern(
        QStringLiteral(R"(\[\s*(%1)\s*,\s*(%1)\s*,\s*(%1)\s*\])")
            .arg(number));

    struct Vertex {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };
    QVector<Vertex> vertices;
    QString verticesText;
    if (extractArrayValue(text, QStringLiteral("vertices"),
                          &verticesText)) {
        auto pointMatches = pointPattern.globalMatch(verticesText);
        while (pointMatches.hasNext()) {
            const QRegularExpressionMatch match = pointMatches.next();
            vertices.push_back({
                match.captured(1).toDouble(),
                match.captured(2).toDouble(),
                match.captured(3).toDouble()
            });
        }
    } else {
        const QRegularExpression verticesSection(
            QStringLiteral(
                R"(^\s*\[vertices\]\s*(.*?)(?=^\s*\[[^\]]+\]|\z))"),
            QRegularExpression::MultilineOption
                | QRegularExpression::DotMatchesEverythingOption);
        const auto sectionMatch = verticesSection.match(text);
        if (sectionMatch.hasMatch()) {
            const QRegularExpression indexedPoint(
                QStringLiteral(
                    R"(^\s*v(\d+)\s*=\s*\[\s*(%1)\s*,\s*(%1)\s*,\s*(%1)\s*\])")
                    .arg(number),
                QRegularExpression::MultilineOption);
            QMap<int, Vertex> indexed;
            auto matches =
                indexedPoint.globalMatch(sectionMatch.captured(1));
            while (matches.hasNext()) {
                const auto match = matches.next();
                indexed.insert(match.captured(1).toInt(),
                               {match.captured(2).toDouble(),
                                match.captured(3).toDouble(),
                                match.captured(4).toDouble()});
            }
            for (const Vertex& vertex : indexed) vertices.push_back(vertex);
        }
    }
    if (vertices.isEmpty()) {
        if (error) {
            *error = QStringLiteral(
                "几何配置中未找到 [vertices] 的 v0、v1... 定义");
        }
        return false;
    }
    if (vertices.size() < 3) {
        if (error) *error = QStringLiteral("vertices 至少需要三个顶点");
        return false;
    }

    const double planeZ = vertices.first().z;
    constexpr double planeTolerance = 1.0e-10;
    GeometrySketch sketch;
    QHash<int, int> vertexToSketch;
    for (int i = 0; i < vertices.size(); ++i) {
        if (std::abs(vertices.at(i).z - planeZ) > planeTolerance) continue;
        vertexToSketch.insert(i, sketch.points.size());
        sketch.points.push_back(
            QPointF(vertices.at(i).x, -vertices.at(i).y));
    }
    if (sketch.points.size() < 3) {
        if (error) *error = QStringLiteral("无法识别几何底面轮廓");
        return false;
    }

    const QRegularExpression edgeSection(
        QStringLiteral(
            R"(^\s*\[edges\]\s*(.*?)(?=^\s*\[[^\]]+\]|\z))"),
        QRegularExpression::MultilineOption
            | QRegularExpression::DotMatchesEverythingOption);
    const QRegularExpression vertexPair(
        QStringLiteral(R"(\bv\s*=\s*\[\s*(\d+)\s*,\s*(\d+)\s*\])"));
    const QRegularExpression middlePoint(
        QStringLiteral(
            R"(\bp\s*=\s*\[\s*(%1)\s*,\s*(%1)\s*,\s*(%1)\s*\])")
            .arg(number));

    auto edgeMatches = edgeSection.globalMatch(text);
    while (edgeMatches.hasNext()) {
        const QString body = edgeMatches.next().captured(1);
        if (!body.contains(
                QRegularExpression(QStringLiteral(R"(\btype\s*=\s*["']?arc)"),
                                   QRegularExpression::CaseInsensitiveOption))) {
            continue;
        }
        const QRegularExpressionMatch pair = vertexPair.match(body);
        const QRegularExpressionMatch middle = middlePoint.match(body);
        if (!pair.hasMatch() || !middle.hasMatch()) continue;
        const int sourceStart = pair.captured(1).toInt();
        const int sourceEnd = pair.captured(2).toInt();
        if (!vertexToSketch.contains(sourceStart)
            || !vertexToSketch.contains(sourceEnd)) {
            continue;
        }
        sketch.arcs.push_back({
            vertexToSketch.value(sourceStart),
            vertexToSketch.value(sourceEnd),
            QPointF(middle.captured(1).toDouble(),
                    -middle.captured(2).toDouble())
        });
    }

    canvas_->setSketch(sketch);
    statusLabel_->setText(
        QStringLiteral("已读取 %1：%2 个底面顶点，%3 条圆弧")
            .arg(info.fileName())
            .arg(sketch.points.size())
            .arg(sketch.arcs.size()));
    return true;
}

} // namespace SF::GUI
