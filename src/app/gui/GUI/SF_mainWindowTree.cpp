/// @file SF_mainWindowTree.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_mainWindow.h"

#include "SF_switchButton.h"
#include "SF_GeometryPanel.h"
#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QSignalBlocker>

#include <QAbstractButton>
#include <QAction>
#include <QBrush>
#include <QColor>
#include <QComboBox>
#include <QFont>
#include <QHash>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QProgressBar>
#include <QRegularExpression>
#include <QToolButton>
#include <QTreeWidget>
#include <QTreeWidgetItemIterator>

#include <algorithm>
#include <cmath>

namespace SF::GUI {
namespace {

QString treeValue(const ConfigEntry& entry)
{
    if (!entry.coordinateTuple) return entry.value;

    QString coordinate = entry.value.trimmed();
    if (coordinate.startsWith(QLatin1Char('('))
        && coordinate.endsWith(QLatin1Char(')'))) {
        coordinate = coordinate.mid(1, coordinate.size() - 2);
    }
    coordinate.replace(QLatin1Char(','), QLatin1Char(' '));
    return coordinate.simplified();
}

QString displayTime(const OutputTaskInfo& task)
{
    if (!task.latestTime.isEmpty()) {
        bool ok = false;
        const double value = task.latestTime.toDouble(&ok);
        if (ok) return QString::number(value, 'e', 2);
        return task.latestTime;
    }
    return QStringLiteral("-");
}

QString displayProgress(double progress)
{
    const int percent = static_cast<int>(
        std::round(std::clamp(progress, 0.0, 1.0) * 100.0));
    return QStringLiteral("%1%").arg(percent);
}

QVector<int> integerList(const QString& text, int expected, int fallback = 0)
{
    QVector<int> values;
    const QRegularExpression number(QStringLiteral(R"([-+]?\d+)"));
    auto match = number.globalMatch(text);
    while (match.hasNext() && values.size() < expected) {
        values.push_back(match.next().captured().toInt());
    }
    while (values.size() < expected) values.push_back(fallback);
    return values;
}

QVector<double> decimalList(const QString& text, int expected, double fallback = 1.0)
{
    QVector<double> values;
    const QRegularExpression number(
        QStringLiteral(R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"));
    auto match = number.globalMatch(text);
    while (match.hasNext() && values.size() < expected) {
        values.push_back(match.next().captured().toDouble());
    }
    while (values.size() < expected) values.push_back(fallback);
    return values;
}

struct BlockTreeValue {
    QString type = QStringLiteral("hex");
    QVector<int> vertices{0, 0, 0, 0, 0, 0, 0, 0};
    QVector<int> seed{1, 1, 1};
    QVector<double> grading{1.0, 1.0, 1.0};
};

BlockTreeValue parseBlockTreeValue(const QString& value)
{
    BlockTreeValue parsed;
    const QRegularExpression pattern(QStringLiteral(
        R"(^\s*(hex|wedge|prism)\s*\(([^)]*)\)\s*\(([^)]*)\)\s*simpleGrading\s*\(([^)]*)\))"),
        QRegularExpression::CaseInsensitiveOption);
    const QRegularExpressionMatch match = pattern.match(value);
    if (!match.hasMatch()) return parsed;
    parsed.type = match.captured(1).toLower();
    parsed.vertices = integerList(match.captured(2), 8);
    parsed.seed = integerList(match.captured(3), 3, 1);
    parsed.grading = decimalList(match.captured(4), 3, 1.0);
    return parsed;
}

QString blockTreeValue(const QString& type, const QVector<int>& vertices,
                       const QVector<double>& seed,
                       const QVector<double>& grading)
{
    return QStringLiteral("%1 (%2 %3 %4 %5 %6 %7 %8 %9) (%10 %11 %12) simpleGrading (%13 %14 %15)")
        .arg(type)
        .arg(vertices.value(0)).arg(vertices.value(1)).arg(vertices.value(2))
        .arg(vertices.value(3)).arg(vertices.value(4)).arg(vertices.value(5))
        .arg(vertices.value(6)).arg(vertices.value(7))
        .arg(seed.value(0)).arg(seed.value(1)).arg(seed.value(2))
        .arg(grading.value(0), 0, 'g', 12).arg(grading.value(1), 0, 'g', 12)
        .arg(grading.value(2), 0, 'g', 12);
}

QString faceTreeValue(const QVector<int>& vertices)
{
    return QStringLiteral("(%1 %2 %3 %4)")
        .arg(vertices.value(0)).arg(vertices.value(1))
        .arg(vertices.value(2)).arg(vertices.value(3));
}

} // namespace

void MainWindow::setWorkspace(const QString& displayName,
                              const QString& absolutePath)
{
    workspaceLabel_->setText(displayName);
    workspaceLabel_->setToolTip(absolutePath);
}

void MainWindow::setProject(const QVector<ConfigDocument>& documents,
                            const QVector<ProjectModule>& modules)
{
    const int previousModule = currentModuleIndex_;
    documents_ = documents;
    modules_ = modules;

    moduleSelectorMenu_->clear();
    for (int i = 0; i < modules_.size(); ++i) {
        const ProjectModule& module = modules_.at(i);
        QAction* action = moduleSelectorMenu_->addAction(
            moduleIcon(module.kind), module.title);
        action->setCheckable(true);
        action->setData(i);
        connect(action, &QAction::triggered, this, [this, i] {
            focusModule(i);
        });
    }
    const int selectedModule = modules_.isEmpty()
        ? -1 : std::clamp(previousModule, 0,
                          static_cast<int>(modules_.size()) - 1);
    currentModuleIndex_ = std::max(0, selectedModule);
    rebuildTree();
    if (selectedModule >= 0) focusModule(selectedModule);
}

void MainWindow::setOutputTasks(const QVector<OutputTaskInfo>& tasks)
{
    outputTasks_ = tasks;
    rebuildTree();
    if (currentModuleIndex_ >= 0 && currentModuleIndex_ < modules_.size()) {
        focusModule(currentModuleIndex_);
    }
}

void MainWindow::rebuildTree()
{
    rebuildingTree_ = true;
    configTree_->clear();
    const int selectedModule = currentModuleIndex_;
    for (int moduleIndex = 0; moduleIndex < modules_.size(); ++moduleIndex) {
        const ProjectModule& module = modules_.at(moduleIndex);
        auto* moduleItem = new QTreeWidgetItem(configTree_);
        moduleItem->setText(0, module.title);
        moduleItem->setIcon(0, moduleIcon(module.kind));
        moduleItem->setData(0, ModuleRole, moduleIndex);
        QFont font = moduleItem->font(0);
        font.setBold(true);
        moduleItem->setFont(0, font);
        const bool containsDisplayedOutput =
            module.kind == ProjectModuleKind::Output
            && std::any_of(outputTasks_.cbegin(), outputTasks_.cend(),
                           [](const OutputTaskInfo& task) {
                               return task.displayed;
                           });
        moduleItem->setExpanded(
            moduleIndex == selectedModule || containsDisplayedOutput);

        if (module.kind == ProjectModuleKind::Output) {
            for (const OutputTaskInfo& task : outputTasks_) {
                auto* taskItem = new QTreeWidgetItem(moduleItem);
                taskItem->setText(0, task.name);
                taskItem->setIcon(0, moduleIcon(ProjectModuleKind::Output));
                taskItem->setData(0, OutputTaskRole, task.name);
                taskItem->setToolTip(
                    0, task.pvdPath.isEmpty()
                           ? QStringLiteral("尚未生成 result/%1.pvd")
                                 .arg(task.name)
                           : task.pvdPath);
                QFont taskFont = taskItem->font(0);
                taskFont.setBold(task.displayed);
                taskItem->setFont(0, taskFont);
                if (task.displayed) {
                    // 这是“右侧正在显示此 PVD”的持久状态，不依赖树的当前选中项。
                    const QBrush displayedBrush(QColor(113, 133, 111, 138));
                    taskItem->setBackground(0, displayedBrush);
                    taskItem->setBackground(1, displayedBrush);
                    taskItem->setForeground(0, QColor(255, 255, 255));
                }
                taskItem->setExpanded(task.displayed);
                auto* progress = new QProgressBar(configTree_);
                progress->setRange(0, 100);
                progress->setValue(static_cast<int>(
                    std::round(std::clamp(task.progress, 0.0, 1.0)
                               * 100.0)));
                progress->setTextVisible(false);
                progress->setFixedHeight(18);
                if (task.displayed) {
                    progress->setStyleSheet(QStringLiteral(
                        "QProgressBar{background:rgba(113,133,111,92);"
                        "border:2px solid #71856f;border-radius:5px;}"
                        "QProgressBar::chunk{background:#71856f;"
                        "border-radius:3px;}"));
                }
                progress->setToolTip(
                    task.displayed
                        ? QStringLiteral("右侧窗口正在显示 %1.pvd · %2帧")
                              .arg(task.name).arg(task.frameCount)
                        : task.frameCount > 0
                        ? QStringLiteral("%1 · %2帧")
                              .arg(displayProgress(task.progress))
                              .arg(task.frameCount)
                        : QStringLiteral("尚未计算"));
                configTree_->setItemWidget(taskItem, 1, progress);

                auto* residualItem = new QTreeWidgetItem(taskItem);
                residualItem->setText(0, QStringLiteral("残差显示"));
                residualItem->setText(1, QStringLiteral("待接入"));
                residualItem->setData(0, OutputTaskRole, task.name);
                residualItem->setData(0, OutputActionRole,
                                      QStringLiteral("residual"));
                residualItem->setForeground(1, QColor(116, 106, 91));

                auto* timeItem = new QTreeWidgetItem(taskItem);
                timeItem->setText(0, QStringLiteral("时间"));
                timeItem->setText(1, displayTime(task));
                timeItem->setForeground(1, QColor(116, 106, 91));

                auto* fileItem = new QTreeWidgetItem(taskItem);
                fileItem->setText(0, QStringLiteral("结果文件"));
                fileItem->setText(1, task.name + QStringLiteral(".pvd"));
                fileItem->setForeground(1, QColor(116, 106, 91));
            }
            continue;
        }

        // 几何模块 section 缓存
        QHash<QString, QTreeWidgetItem*> geomSections_;
        auto ensureGeomSection = [&](const QString& name) -> QTreeWidgetItem* {
            if (!geomSections_.contains(name)) {
                auto* sec = new QTreeWidgetItem(moduleItem);
                sec->setText(0, name);
                QFont sf = sec->font(0); sf.setBold(true); sec->setFont(0, sf);
                sec->setExpanded(true);
                geomSections_[name] = sec;
            }
            return geomSections_[name];
        };

        // 块缓冲：按 sectionLabel 分组 [blocks] 条目（用于后续批量渲染）
        QHash<QString, QVector<ProjectEntryRef>> blockBuf_;
        // Patch 缓冲：按 section 名分组（boundary/* 条目）
        QHash<QString, QVector<ProjectEntryRef>> patchBuf_;
        QStringList patchNames;

        QHash<QString, QStringList> sectionLabels;
        int blockMeshDocumentIndex = -1;
        for (const ProjectEntryRef& ref : module.entries) {
            if (ref.documentIndex < 0 || ref.documentIndex >= documents_.size()) {
                continue;
            }
            if (documents_.at(ref.documentIndex).fileName().compare(
                    QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0) {
                blockMeshDocumentIndex = ref.documentIndex;
            }
            const auto& entries = documents_.at(ref.documentIndex).entries();
            if (ref.entryIndex < 0 || ref.entryIndex >= entries.size()) continue;
            const ConfigEntry& entry = entries.at(ref.entryIndex);
            if (entry.caseSection) continue;
            const QString sectionKey =
                QStringLiteral("%1:%2")
                    .arg(ref.documentIndex)
                    .arg(entry.section.toLower());
            if (!sectionLabels[sectionKey].contains(entry.sectionLabel)) {
                sectionLabels[sectionKey].push_back(entry.sectionLabel);
            }
        }

        // 即使集合暂时为空也保留三组根节点，删除最后一个条目后仍可右键继续新增。
        if (blockMeshDocumentIndex >= 0) {
            const auto initGeometrySection = [&](const QString& title,
                                                 const QString& sectionName) {
                auto* section = ensureGeomSection(title);
                section->setData(0, DocumentRole, blockMeshDocumentIndex);
                section->setData(0, SectionRole, sectionName);
            };
            initGeometrySection(QStringLiteral("Vertices"), QStringLiteral("vertices"));
            initGeometrySection(QStringLiteral("Blocks"), QStringLiteral("blocks"));
            initGeometrySection(QStringLiteral("Patches"), QStringLiteral("boundary"));
        }

        QHash<QString, QTreeWidgetItem*> sectionItems;
        QHash<QString, QTreeWidgetItem*> occurrenceItems;
        QHash<QString, QTreeWidgetItem*> currentRecordItems;
        QHash<QString, QTreeWidgetItem*> boundaryTargetItems;
        QHash<QString, QTreeWidgetItem*> boundaryFieldItems;
        if (module.kind == ProjectModuleKind::BoundaryCondition) {
            for (const QString& target : module.targetNames) {
                auto* targetItem = new QTreeWidgetItem(moduleItem);
                targetItem->setText(0, target);
                targetItem->setData(0, SectionRole,
                                    QStringLiteral("boundaryTarget/%1").arg(target));
                QFont targetFont = targetItem->font(0);
                targetFont.setWeight(QFont::DemiBold);
                targetItem->setFont(0, targetFont);
                boundaryTargetItems.insert(target.toLower(), targetItem);
            }
        }

        QString turbulenceType;
        if (module.kind == ProjectModuleKind::Turbulence) {
            for (const ProjectEntryRef& ref : module.entries) {
                if (ref.documentIndex < 0 || ref.documentIndex >= documents_.size()) continue;
                const auto& entries = documents_.at(ref.documentIndex).entries();
                if (ref.entryIndex < 0 || ref.entryIndex >= entries.size()) continue;
                const ConfigEntry& candidate = entries.at(ref.entryIndex);
                if (candidate.key.compare(QStringLiteral("simulationType"),
                                          Qt::CaseInsensitive) == 0) {
                    turbulenceType = candidate.value.trimmed();
                    break;
                }
            }
        }

        auto configureEntryItem =
            [this](QTreeWidgetItem* entryItem,
                   const ProjectEntryRef& ref,
                   const ConfigEntry& entry) {
                entryItem->setData(0, DocumentRole, ref.documentIndex);
                entryItem->setData(0, EntryRole, ref.entryIndex);
                entryItem->setToolTip(0, entry.value);
                const QStringList choices = enumChoices(entry);
                if (!choices.isEmpty()) {
                    auto* selector = new QComboBox(configTree_);
                    selector->addItems(choices);
                    int current = choices.indexOf(entry.value.trimmed());
                    if (current < 0) {
                        for (int i = 0; i < choices.size(); ++i) {
                            if (choices.at(i).compare(entry.value.trimmed(),
                                                      Qt::CaseInsensitive) == 0) {
                                current = i;
                                break;
                            }
                        }
                    }
                    selector->setCurrentIndex(std::max(0, current));
                    configTree_->setItemWidget(entryItem, 1, selector);
                    connect(selector, &QComboBox::currentTextChanged, this,
                            [this, ref](const QString& value) {
                        emit configValueEdited(ref.documentIndex, ref.entryIndex,
                                               value);
                    });
                } else if (entry.isBoolean()) {
                    auto* toggle = new SwitchButton(configTree_);
                    toggle->setChecked(entry.value.trimmed().compare(
                                           QStringLiteral("true"),
                                           Qt::CaseInsensitive) == 0);
                    configTree_->setItemWidget(entryItem, 1, toggle);
                    connect(toggle, &QAbstractButton::toggled, this,
                            [this, ref](bool checked) {
                        emit configValueEdited(
                            ref.documentIndex, ref.entryIndex,
                            checked ? QStringLiteral("true")
                                    : QStringLiteral("false"));
                    });
                } else {
                    const QString value = treeValue(entry);
                    entryItem->setText(1, value);
                    entryItem->setData(1, OriginalValueRole, value);
                    entryItem->setToolTip(1, entry.value);
                    entryItem->setForeground(1, QColor(116, 106, 91));
                    entryItem->setFlags(entryItem->flags()
                                        | Qt::ItemIsEditable);
                }
            };

        for (const ProjectEntryRef& ref : module.entries) {
            if (ref.documentIndex < 0 || ref.documentIndex >= documents_.size()) {
                continue;
            }
            const auto& entries = documents_.at(ref.documentIndex).entries();
            if (ref.entryIndex < 0 || ref.entryIndex >= entries.size()) continue;
            const ConfigEntry& entry = entries.at(ref.entryIndex);

            if (module.kind == ProjectModuleKind::Turbulence
                && entry.section.startsWith(QStringLiteral("RAS"),
                                            Qt::CaseInsensitive)
                && turbulenceType.compare(QStringLiteral("RAS"),
                                          Qt::CaseInsensitive) != 0) {
                continue;
            }

            if (module.kind == ProjectModuleKind::InitialCondition
                && entry.key.compare(QStringLiteral("internalField"),
                                     Qt::CaseInsensitive) == 0) {
                auto* fieldItem = new QTreeWidgetItem(moduleItem);
                fieldItem->setText(0, documents_.at(ref.documentIndex).fileName());
                configureEntryItem(fieldItem, ref, entry);
                continue;
            }

            if (module.kind == ProjectModuleKind::BoundaryCondition
                && entry.section.startsWith(QStringLiteral("boundaryField/"),
                                            Qt::CaseInsensitive)) {
                const QString patchName = entry.section.section(QLatin1Char('/'), 1, 1);
                QTreeWidgetItem* patchItem = boundaryTargetItems.value(patchName.toLower());
                if (!patchItem) {
                    patchItem = new QTreeWidgetItem(moduleItem);
                    patchItem->setText(0, patchName);
                    boundaryTargetItems.insert(patchName.toLower(), patchItem);
                }
                const QString fieldName = documents_.at(ref.documentIndex).fileName();
                const QString fieldKey = patchName.toLower() + QLatin1Char('/')
                    + fieldName.toLower();
                QTreeWidgetItem* fieldItem = boundaryFieldItems.value(fieldKey);
                if (!fieldItem) {
                    fieldItem = new QTreeWidgetItem(patchItem);
                    fieldItem->setText(0, fieldName);
                    boundaryFieldItems.insert(fieldKey, fieldItem);
                }
                auto* entryItem = new QTreeWidgetItem(fieldItem);
                entryItem->setText(0, entry.key);
                configureEntryItem(entryItem, ref, entry);
                continue;
            }

            if (entry.caseSection) {
                const QString csLower = entry.section.toLower();
                // === vertices: ParaView风格三坐标输入 ===
                if (csLower == QStringLiteral("vertices")) {
                    auto* sec = ensureGeomSection(QStringLiteral("Vertices"));
                    sec->setData(0, DocumentRole, ref.documentIndex);
                    sec->setData(0, SectionRole, QStringLiteral("vertices"));
                    auto* vItem = new QTreeWidgetItem(sec);
                    vItem->setText(0, QString());
                    vItem->setData(0, DocumentRole, ref.documentIndex);
                    vItem->setData(0, EntryRole, ref.entryIndex);
                    vItem->setData(0, SectionRole, QStringLiteral("vertices"));
                    vItem->setFirstColumnSpanned(true);
                    double px=0,py=0,pz=0;
                    QString val = entry.value.trimmed();
                    val.remove(QLatin1Char('('));val.remove(QLatin1Char(')'));
                    val.remove(QLatin1Char('['));val.remove(QLatin1Char(']'));
                    val.replace(QLatin1Char(','),QLatin1Char(' '));
                    const QStringList sp = val.simplified().split(QLatin1Char(' '));
                    if (sp.size()>=3){px=sp[0].toDouble();py=sp[1].toDouble();pz=sp[2].toDouble();}
                    auto* w = new QWidget(configTree_);
                    auto* hr = new QHBoxLayout(w);hr->setContentsMargins(0,2,0,2);hr->setSpacing(3);
                    auto* label = new QLabel(entry.key, w);
                    label->setFixedWidth(30);
                    label->setStyleSheet(QStringLiteral("color:#4a443a;font-weight:600;"));
                    auto mk=[w](double v)->QDoubleSpinBox*{
                        auto* s=new QDoubleSpinBox(w);s->setRange(-1e9,1e9);s->setDecimals(6);
                        s->setFixedWidth(58);s->setButtonSymbols(QAbstractSpinBox::NoButtons);s->setValue(v);
                        s->setStyleSheet("QDoubleSpinBox{background:#fffdf8;border:1px solid #aa9f8f;border-radius:4px;padding:1px 3px;color:#3d382f;font-size:11px;}");
                        return s;};
                    auto*sx=mk(px);auto*sy=mk(py);auto*sz=mk(pz);
                    hr->addWidget(label);hr->addWidget(sx);hr->addWidget(sy);hr->addWidget(sz);hr->addStretch();
                    configTree_->setItemWidget(vItem,0,w);
                    auto save=[this,ref,sx,sy,sz](){emit configValueEdited(ref.documentIndex,ref.entryIndex,QStringLiteral("[%1,%2,%3]").arg(sx->value(),0,'f',6).arg(sy->value(),0,'f',6).arg(sz->value(),0,'f',6));};
                    connect(sx,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
                    connect(sy,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
                    connect(sz,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
                    continue;
                }
                // === blocks: 缓冲后统一渲染 ===
                if (csLower == QStringLiteral("blocks")) {
                    const QString occKey = entry.sectionLabel;
                    blockBuf_[occKey].push_back(ref);
                    continue;
                }
                // === edges ===
                if (csLower == QStringLiteral("edges")) {
                    auto* sec = ensureGeomSection(QStringLiteral("Edges"));
                    auto* entryItem = new QTreeWidgetItem(sec);
                    entryItem->setText(0, entryLabel(entry));
                    configureEntryItem(entryItem, ref, entry);
                    continue;
                }
                // === boundary patch: 缓冲后统一渲染 ===
                if (csLower.startsWith(QStringLiteral("boundary/"))) {
                    patchBuf_[entry.section].push_back(ref);
                    if (!patchNames.contains(entry.section)) patchNames.push_back(entry.section);
                    continue;
                }
                // 默认
                auto* entryItem = new QTreeWidgetItem(moduleItem);
                entryItem->setText(0, entryLabel(entry));
                configureEntryItem(entryItem, ref, entry);
                continue;
            }

            const QString sectionKey =
                QStringLiteral("%1:%2")
                    .arg(ref.documentIndex)
                    .arg(entry.section.toLower());
            QTreeWidgetItem* sectionItem = entry.section.isEmpty()
                ? moduleItem : sectionItems.value(sectionKey);
            if (!sectionItem) {
                sectionItem = new QTreeWidgetItem(moduleItem);
                sectionItem->setText(
                    0, QStringLiteral("[%1]").arg(entry.section));
                sectionItem->setData(0, DocumentRole, ref.documentIndex);
                sectionItem->setData(0, SectionRole, entry.section);
                QFont sectionFont = sectionItem->font(0);
                sectionFont.setWeight(QFont::DemiBold);
                sectionItem->setFont(0, sectionFont);
                sectionItems.insert(sectionKey, sectionItem);
            }

            QTreeWidgetItem* parentItem = sectionItem;
            if (sectionLabels.value(sectionKey).size() > 1) {
                const QString occurrenceKey =
                    sectionKey + QLatin1Char(':') + entry.sectionLabel;
                parentItem = occurrenceItems.value(occurrenceKey);
                if (!parentItem) {
                    parentItem = new QTreeWidgetItem(sectionItem);
                    const int occurrence =
                        sectionLabels.value(sectionKey)
                            .indexOf(entry.sectionLabel) + 1;
                    parentItem->setText(
                        0, QStringLiteral("%1 %2")
                               .arg(entry.section)
                               .arg(occurrence));
                    occurrenceItems.insert(occurrenceKey, parentItem);
                }
            } else if (entry.key.compare(QStringLiteral("set"),
                                         Qt::CaseInsensitive) == 0
                       || entry.key.startsWith(QStringLiteral("set"),
                                               Qt::CaseInsensitive)) {
                QString recordName = entry.value.trimmed();
                if ((recordName.startsWith(QLatin1Char('"'))
                     && recordName.endsWith(QLatin1Char('"')))
                    || (recordName.startsWith(QLatin1Char('\''))
                        && recordName.endsWith(QLatin1Char('\'')))) {
                    recordName = recordName.mid(1, recordName.size() - 2);
                }
                auto* recordItem = new QTreeWidgetItem(sectionItem);
                recordItem->setText(0, recordName);
                configureEntryItem(recordItem, ref, entry);
                recordItem->setText(1, QString());
                currentRecordItems.insert(sectionKey, recordItem);
                continue;
            } else if (currentRecordItems.contains(sectionKey)) {
                parentItem = currentRecordItems.value(sectionKey);
            }

            // 普通 key=value 节点（含 blockMeshDict 格式几何条目）
            const QString secLower2 = entry.section.toLower();
            // === vertices（任何格式） ===
            if (secLower2 == QStringLiteral("vertices")) {
                auto* vItem = new QTreeWidgetItem(parentItem);
                vItem->setText(0, QString());
                vItem->setData(0, DocumentRole, ref.documentIndex);
                vItem->setData(0, EntryRole, ref.entryIndex);
                vItem->setData(0, SectionRole, QStringLiteral("vertices"));
                vItem->setFirstColumnSpanned(true);
                double px=0,py=0,pz=0;
                QString c = entry.value.trimmed();
                c.remove(QLatin1Char('('));c.remove(QLatin1Char(')'));
                c.remove(QLatin1Char('['));c.remove(QLatin1Char(']'));
                c.replace(QLatin1Char(','),QLatin1Char(' '));
                QStringList parts = c.simplified().split(QLatin1Char(' '));
                if (parts.size()>=3){px=parts[0].toDouble();py=parts[1].toDouble();pz=parts[2].toDouble();}
                auto* w = new QWidget(configTree_);
                auto* hr = new QHBoxLayout(w);
                hr->setContentsMargins(0,2,0,2);hr->setSpacing(3);
                auto* label = new QLabel(entryLabel(entry), w);
                label->setFixedWidth(30);
                label->setStyleSheet(QStringLiteral("color:#4a443a;font-weight:600;"));
                auto mk=[w](double v)->QDoubleSpinBox*{
                    auto* s=new QDoubleSpinBox(w);s->setRange(-1e9,1e9);s->setDecimals(6);
                    s->setFixedWidth(58);s->setButtonSymbols(QAbstractSpinBox::NoButtons);s->setValue(v);
                    s->setStyleSheet("QDoubleSpinBox{background:#fffdf8;border:1px solid #aa9f8f;border-radius:4px;padding:1px 3px;color:#3d382f;font-size:11px;}");
                    return s;};
                auto*sx=mk(px);auto*sy=mk(py);auto*sz=mk(pz);
                hr->addWidget(label);hr->addWidget(sx);hr->addWidget(sy);hr->addWidget(sz);hr->addStretch();
                configTree_->setItemWidget(vItem,0,w);
                auto save=[this,ref,sx,sy,sz](){emit configValueEdited(ref.documentIndex,ref.entryIndex,QStringLiteral("[%1,%2,%3]").arg(sx->value(),0,'f',6).arg(sy->value(),0,'f',6).arg(sz->value(),0,'f',6));};
                connect(sx,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
                connect(sy,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
                connect(sz,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,save);
            // === blocks（blockMeshDict格式：缓冲后统一渲染）===
            } else if (secLower2 == QStringLiteral("blocks")) {
                blockBuf_[entry.sectionLabel].push_back(ref);
            // === edges ===
            } else if (secLower2 == QStringLiteral("edges")) {
                auto* entryItem = new QTreeWidgetItem(parentItem);
                entryItem->setText(0, entryLabel(entry));
                configureEntryItem(entryItem, ref, entry);
            // === boundary / boundary/*（blockMeshDict 或 .sf 格式）===
            } else if (secLower2 == QStringLiteral("boundary")
                       || secLower2.startsWith(QStringLiteral("boundary/"))) {
                patchBuf_[entry.section].push_back(ref);
                if (!patchNames.contains(entry.section)) patchNames.push_back(entry.section);
            } else {
                auto* entryItem = new QTreeWidgetItem(parentItem);
                entryItem->setText(0, entryLabel(entry));
                configureEntryItem(entryItem, ref, entry);
            }
        }

        // blockMeshDict 的 block 行在此一次性展开为形状、两行顶点、seed、grading。
        for (auto it = blockBuf_.cbegin(); it != blockBuf_.cend(); ++it) {
            if (it.value().isEmpty()) continue;
            const int documentIndex = it.value().first().documentIndex;
            int pointCount = 0;
            if (documentIndex >= 0 && documentIndex < documents_.size()) {
                for (const ConfigEntry& candidate : documents_.at(documentIndex).entries()) {
                    if (candidate.section.compare(QStringLiteral("vertices"),
                                                  Qt::CaseInsensitive) == 0) {
                        ++pointCount;
                    }
                }
            }
            auto* section = ensureGeomSection(QStringLiteral("Blocks"));
            section->setData(0, DocumentRole, documentIndex);
            section->setData(0, SectionRole, QStringLiteral("blocks"));
            for (const ProjectEntryRef& ref : it.value()) {
                if (ref.documentIndex < 0 || ref.documentIndex >= documents_.size()) continue;
                const auto& entries = documents_.at(ref.documentIndex).entries();
                if (ref.entryIndex < 0 || ref.entryIndex >= entries.size()) continue;
                const ConfigEntry& entry = entries.at(ref.entryIndex);
                const BlockTreeValue parsed = parseBlockTreeValue(entry.value);

                auto* blockItem = new QTreeWidgetItem(section);
                blockItem->setText(0, entry.key);
                blockItem->setData(0, DocumentRole, ref.documentIndex);
                blockItem->setData(0, EntryRole, ref.entryIndex);
                blockItem->setData(0, SectionRole, QStringLiteral("blocks"));
                int blockId = -1;
                for (int index = 0; index <= ref.entryIndex; ++index) {
                    if (entries.at(index).section.compare(QStringLiteral("blocks"),
                                                          Qt::CaseInsensitive) == 0) {
                        ++blockId;
                    }
                }
                blockItem->setData(0, BlockIDRole, blockId);
                blockItem->setExpanded(true);
                auto* typeSelector = new CircleRadioGroup(
                    {QStringLiteral("hex"), QStringLiteral("wedge"),
                     QStringLiteral("prism")}, configTree_);
                typeSelector->setCurrentIndex(
                    parsed.type == QStringLiteral("wedge") ? 1
                    : parsed.type == QStringLiteral("prism") ? 2 : 0);
                configTree_->setItemWidget(blockItem, 1, typeSelector);

                auto* lowerVertices = new QTreeWidgetItem(section);
                auto* upperVertices = new QTreeWidgetItem(section);
                lowerVertices->setFirstColumnSpanned(true);
                upperVertices->setFirstColumnSpanned(true);
                auto* lowerRow = new VerticesRow(pointCount, configTree_);
                auto* upperRow = new VerticesRow(pointCount, configTree_);
                lowerRow->setValues(parsed.vertices[0], parsed.vertices[1],
                                    parsed.vertices[2], parsed.vertices[3]);
                upperRow->setValues(parsed.vertices[4], parsed.vertices[5],
                                    parsed.vertices[6], parsed.vertices[7]);
                configTree_->setItemWidget(lowerVertices, 0, lowerRow);
                configTree_->setItemWidget(upperVertices, 0, upperRow);

                auto* seedItem = new QTreeWidgetItem(section);
                auto* gradingItem = new QTreeWidgetItem(section);
                seedItem->setFirstColumnSpanned(true);
                gradingItem->setFirstColumnSpanned(true);
                auto* seedRow = new TripleSpinRow(QStringLiteral("seed:"), true,
                                                   configTree_);
                auto* gradingRow = new TripleSpinRow(QStringLiteral("grading:"), false,
                                                      configTree_);
                seedRow->setValues(parsed.seed[0], parsed.seed[1], parsed.seed[2]);
                gradingRow->setValues(parsed.grading[0], parsed.grading[1],
                                      parsed.grading[2]);
                configTree_->setItemWidget(seedItem, 0, seedRow);
                configTree_->setItemWidget(gradingItem, 0, gradingRow);

                auto saveBlock = [this, ref, typeSelector, lowerRow, upperRow,
                                  seedRow, gradingRow] {
                    QVector<int> vertices = lowerRow->values();
                    const QVector<int> top = upperRow->values();
                    vertices += top;
                    emit configValueEdited(
                        ref.documentIndex, ref.entryIndex,
                        blockTreeValue(typeSelector->currentIndex() == 1
                                           ? QStringLiteral("wedge")
                                           : typeSelector->currentIndex() == 2
                                               ? QStringLiteral("prism")
                                               : QStringLiteral("hex"),
                                       vertices, seedRow->values(), gradingRow->values()));
                };
                connect(typeSelector, &CircleRadioGroup::currentIndexChanged,
                        this, saveBlock);
                connect(lowerRow, &VerticesRow::valuesChanged, this, saveBlock);
                connect(upperRow, &VerticesRow::valuesChanged, this, saveBlock);
                connect(seedRow, &TripleSpinRow::valuesChanged, this, saveBlock);
                connect(gradingRow, &TripleSpinRow::valuesChanged, this, saveBlock);
            }
        }

        // blockMeshDict 的 boundary 被呈现为 Patches 集合；每个面由四个顶点下拉框组成。
        if (!patchNames.isEmpty()) {
            auto* section = ensureGeomSection(QStringLiteral("Patches"));
            const int documentIndex = patchBuf_.value(patchNames.first()).first().documentIndex;
            section->setData(0, DocumentRole, documentIndex);
            section->setData(0, SectionRole, QStringLiteral("boundary"));

            int pointCount = 0;
            if (documentIndex >= 0 && documentIndex < documents_.size()) {
                for (const ConfigEntry& candidate : documents_.at(documentIndex).entries()) {
                    if (candidate.section.compare(QStringLiteral("vertices"),
                                                  Qt::CaseInsensitive) == 0) ++pointCount;
                }
            }
            for (int patchId = 0; patchId < patchNames.size(); ++patchId) {
                const QString& patchSection = patchNames.at(patchId);
                const QVector<ProjectEntryRef>& refs = patchBuf_.value(patchSection);
                if (refs.isEmpty()) continue;
                const QString patchName = patchSection.section(QLatin1Char('/'), -1);
                ProjectEntryRef typeRef;
                QVector<ProjectEntryRef> faceRefs;
                for (const ProjectEntryRef& ref : refs) {
                    const ConfigEntry& entry = documents_.at(ref.documentIndex).entries().at(ref.entryIndex);
                    if (entry.key.compare(QStringLiteral("type"), Qt::CaseInsensitive) == 0)
                        typeRef = ref;
                    else if (entry.key.startsWith(QStringLiteral("face"), Qt::CaseInsensitive))
                        faceRefs.push_back(ref);
                }
                auto* patchItem = new QTreeWidgetItem(section);
                patchItem->setData(0, DocumentRole, documentIndex);
                patchItem->setData(0, SectionRole, QStringLiteral("boundary"));
                patchItem->setData(0, PatchNameRole, patchName);
                patchItem->setData(0, PatchIDRole, patchId);
                auto* nameEditor = new QLineEdit(patchName, configTree_);
                nameEditor->setFrame(false);
                nameEditor->setStyleSheet(QStringLiteral("QLineEdit{color:#4a443a;font-weight:600;}"));
                configTree_->setItemWidget(patchItem, 0, nameEditor);
                auto* controls = new QWidget(configTree_);
                auto* controlsLayout = new QHBoxLayout(controls);
                controlsLayout->setContentsMargins(0, 0, 0, 0);
                controlsLayout->setSpacing(3);
                auto* patchType = new CircleRadioGroup(
                    {QStringLiteral("patch"), QStringLiteral("wall")}, controls);
                bool isWall = false;
                if (typeRef.documentIndex >= 0) {
                    const QString value = documents_.at(typeRef.documentIndex)
                                              .entries().at(typeRef.entryIndex).value;
                    isWall = value.contains(QRegularExpression(
                        QStringLiteral(R"(\bwall\b)"),
                        QRegularExpression::CaseInsensitiveOption));
                }
                patchType->setCurrentIndex(isWall ? 1 : 0);
                controlsLayout->addWidget(patchType, 1);
                configTree_->setItemWidget(patchItem, 1, controls);
                connect(nameEditor, &QLineEdit::cursorPositionChanged, this,
                        [this, patchItem](int, int) {
                    configTree_->setCurrentItem(patchItem);
                    patchItem->setSelected(true);
                });
                connect(nameEditor, &QLineEdit::editingFinished, this,
                        [this, documentIndex, patchName, nameEditor] {
                    emit blockMeshPatchRenameRequested(documentIndex, patchName,
                                                       nameEditor->text().trimmed());
                });
                if (typeRef.documentIndex >= 0) {
                    connect(patchType, &CircleRadioGroup::currentIndexChanged, this,
                            [this, typeRef](int index) {
                        emit configValueEdited(typeRef.documentIndex, typeRef.entryIndex,
                                               index == 0 ? QStringLiteral("type patch;")
                                                          : QStringLiteral("type wall;"));
                    });
                }
                for (int faceId = 0; faceId < faceRefs.size(); ++faceId) {
                    const ProjectEntryRef& faceRef = faceRefs.at(faceId);
                    const ConfigEntry& entry = documents_.at(faceRef.documentIndex)
                                                   .entries().at(faceRef.entryIndex);
                    const QVector<int> vertices = integerList(entry.value, 4);
                    // 面输入与 patch 标题同级，避免树的最后一级再向右缩进。
                    auto* faceItem = new QTreeWidgetItem(section);
                    faceItem->setData(0, DocumentRole, faceRef.documentIndex);
                    faceItem->setData(0, EntryRole, faceRef.entryIndex);
                    faceItem->setData(0, SectionRole, QStringLiteral("boundaryFace"));
                    faceItem->setData(0, PatchIDRole, patchId);
                    faceItem->setData(0, FaceIDRole, faceId);
                    faceItem->setFirstColumnSpanned(true);
                    // Patch 的面只显示四个顶点选择器，不额外占用 face1/face2 标签空间。
                    faceItem->setText(0, QString());
                    auto* faceControls = new QWidget(configTree_);
                    auto* faceLayout = new QHBoxLayout(faceControls);
                    faceLayout->setContentsMargins(0, 0, 0, 0);
                    faceLayout->setSpacing(3);
                    auto* row = new VerticesRow(pointCount, faceControls);
                    row->setValues(vertices[0], vertices[1], vertices[2], vertices[3]);
                    faceLayout->addWidget(row, 1);
                    configTree_->setItemWidget(faceItem, 0, faceControls);
                    connect(row, &VerticesRow::interactionStarted, this,
                            [this, faceItem] {
                        configTree_->setCurrentItem(faceItem);
                        faceItem->setSelected(true);
                    });
                    connect(row, &VerticesRow::valuesChanged, this, [this, faceRef, row] {
                        emit configValueEdited(faceRef.documentIndex, faceRef.entryIndex,
                                               faceTreeValue(row->values()));
                    });
                }
            }
        }
    }
    updateEnabledMenu();
    rebuildingTree_ = false;
}

QString MainWindow::entryLabel(const ConfigEntry& entry) const
{
    if (!entry.caseSection) return entry.key;
    if (entry.section.compare(QStringLiteral("vertices"),
                              Qt::CaseInsensitive) == 0
        || entry.section.compare(QStringLiteral("blocks"),
                                 Qt::CaseInsensitive) == 0
        || entry.section.compare(QStringLiteral("edges"),
                                 Qt::CaseInsensitive) == 0
        || entry.section.startsWith(QStringLiteral("boundary/"),
                                    Qt::CaseInsensitive)) {
        return entry.key;
    }
    if (entry.section.compare(QStringLiteral("mesh"),
                              Qt::CaseInsensitive) == 0) {
        return QStringLiteral("网格文件");
    }
    if (entry.section.compare(QStringLiteral("job"),
                              Qt::CaseInsensitive) == 0) {
        return QStringLiteral("文件名称");
    }
    if (entry.section.compare(QStringLiteral("outputDirectory"),
                              Qt::CaseInsensitive) == 0) {
        return QStringLiteral("输出目录");
    }
    return entry.section;
}



void MainWindow::editTreeItem(QTreeWidgetItem* item)
{
    if (!item) return;
    if (handleOutputTreeClick(item)) return;
    if (!item->data(0, DocumentRole).isValid()) return;
    const int documentIndex = item->data(0, DocumentRole).toInt();
    const int entryIndex = item->data(0, EntryRole).toInt();
    if (documentIndex < 0 || documentIndex >= documents_.size()) return;
    const auto& entries = documents_.at(documentIndex).entries();
    if (entryIndex < 0 || entryIndex >= entries.size()) return;
    if (entries.at(entryIndex).isBoolean()) return;
    configTree_->setCurrentItem(item, 1);
    configTree_->editItem(item, 1);
}

void MainWindow::commitTreeItemEdit(QTreeWidgetItem* item, int column)
{
    if (rebuildingTree_ || !item || column != 1) return;
    if (!item->data(0, DocumentRole).isValid()
        || !item->data(0, EntryRole).isValid()) return;
    const int documentIndex = item->data(0, DocumentRole).toInt();
    const int entryIndex = item->data(0, EntryRole).toInt();
    if (documentIndex < 0 || documentIndex >= documents_.size()) return;
    const auto& entries = documents_.at(documentIndex).entries();
    if (entryIndex < 0 || entryIndex >= entries.size()
        || entries.at(entryIndex).isBoolean()) return;
    const QString original = item->data(1, OriginalValueRole).toString();
    const QString edited = item->text(1).trimmed();
    if (edited == original) return;
    emit configValueEdited(documentIndex, entryIndex, edited);
}

bool MainWindow::handleOutputTreeClick(QTreeWidgetItem* item)
{
    if (!item) return false;
    const QString action = item->data(0, OutputActionRole).toString();
    const QString taskName = item->data(0, OutputTaskRole).toString();
    if (action == QStringLiteral("new")) { requestCreateOutputTask(); return true; }
    if (action == QStringLiteral("residual") && !taskName.isEmpty()) {
        emit outputResidualRequested(taskName); return true;
    }
    if (!taskName.isEmpty()) { emit outputTaskSelected(taskName); return true; }
    return false;
}

void MainWindow::requestCreateOutputTask()
{
    bool accepted = false;
    const QString taskName = QInputDialog::getText(
        this, QStringLiteral("新建输出任务"), QStringLiteral("任务名称"),
        QLineEdit::Normal, QStringLiteral("Sod"), &accepted);
    if (!accepted || taskName.trimmed().isEmpty()) return;
    emit outputTaskCreateRequested(taskName.trimmed());
}

QString MainWindow::currentOutputTaskName() const
{
    for (const OutputTaskInfo& task : outputTasks_)
        if (task.active) return task.name;
    return outputTasks_.isEmpty() ? QString() : outputTasks_.first().name;
}

QStringList MainWindow::enumChoices(const ConfigEntry& entry) const
{
    if (entry.key.compare(QStringLiteral("simulationType"),
                          Qt::CaseInsensitive) == 0) {
        return {QStringLiteral("laminar"), QStringLiteral("DNS"),
                QStringLiteral("LES"), QStringLiteral("RAS")};
    }
    if (entry.key.compare(QStringLiteral("RASModel"),
                          Qt::CaseInsensitive) == 0) {
        return {QStringLiteral("kOmegaSST"), QStringLiteral("kEpsilon"),
                QStringLiteral("RNGkEpsilon"), QStringLiteral("realizableKE"),
                QStringLiteral("SpalartAllmaras")};
    }
    if (entry.key.compare(QStringLiteral("type"), Qt::CaseInsensitive) != 0) return {};
    const QString section = entry.section.toLower();
    if (section == QStringLiteral("velocity") || section == QStringLiteral("p")
        || section == QStringLiteral("rho") || section == QStringLiteral("k")
        || section == QStringLiteral("epsilon") || section == QStringLiteral("omega"))
        return {QStringLiteral("FIXED_VALUE"), QStringLiteral("ZERO_GRADIENT"),
                QStringLiteral("SYMMETRY"), QStringLiteral("EMPTY")};
    if (section == QStringLiteral("output"))
        return {QStringLiteral("time"), QStringLiteral("step")};
    return {};
}

void MainWindow::focusModule(int moduleIndex)
{
    if (moduleIndex < 0 || moduleIndex >= modules_.size()) return;
    currentModuleIndex_ = moduleIndex;
    syncModuleSelector(moduleIndex);
    for (int i = 0; i < configTree_->topLevelItemCount(); ++i) {
        auto* item = configTree_->topLevelItem(i);
        if (item->data(0, ModuleRole).toInt() == moduleIndex) {
            configTree_->scrollToItem(item);
            configTree_->setCurrentItem(item);
            break;
        }
    }
}

void MainWindow::focusMeshModule()
{
    focusModule(ProjectModuleKind::Mesh);
    showVisualizationView();
}

void MainWindow::focusModule(ProjectModuleKind kind)
{
    for (int i = 0; i < modules_.size(); ++i) {
        if (modules_.at(i).kind == kind) { focusModule(i); return; }
    }
}

void MainWindow::syncModuleSelector(int moduleIndex)
{
    if (moduleIndex < 0 || moduleIndex >= modules_.size()) return;
    currentModuleIndex_ = moduleIndex;
    const ProjectModule& module = modules_.at(moduleIndex);
    moduleSelector_->setText(module.title);
    moduleSelector_->setIcon(moduleIcon(module.kind));
}

void MainWindow::updateEnabledMenu()
{
    enabledMenu_->clear();
    for (int di = 0; di < documents_.size(); ++di) {
        const auto& entries = documents_.at(di).entries();
        for (int ei = 0; ei < entries.size(); ++ei) {
            if (entries.at(ei).isBoolean()
                && entries.at(ei).value.trimmed().compare(
                    QStringLiteral("true"), Qt::CaseInsensitive) == 0) {
                QAction* action = enabledMenu_->addAction(
                    QStringLiteral("[%1] %2")
                        .arg(entries.at(ei).section, entries.at(ei).key));
                connect(action, &QAction::triggered, this,
                        [this, di, ei] { selectEntry(di, ei); });
            }
        }
    }
    enabledMenu_->setEnabled(!enabledMenu_->isEmpty());
}

void MainWindow::selectEntry(int documentIndex, int entryIndex)
{
    for (int i = 0; i < configTree_->topLevelItemCount(); ++i) {
        auto* item = configTree_->topLevelItem(i);
        for (int j = 0; j < item->childCount(); ++j) {
            auto* child = item->child(j);
            if (child->data(0, DocumentRole).toInt() == documentIndex
                && child->data(0, EntryRole).toInt() == entryIndex) {
                configTree_->scrollToItem(child);
                configTree_->setCurrentItem(child);
                return;
            }
        }
    }
}

void MainWindow::showModuleContextMenu(const QPoint& position)
{
    QTreeWidgetItem* sectionItem = configTree_->itemAt(position);
    QMenu menu;
    int documentIndex = -1;
    QString section;
    const ProjectModule* module = nullptr;

    if (sectionItem) {
        documentIndex = sectionItem->data(0, DocumentRole).toInt();
        section = sectionItem->data(0, SectionRole).toString();
        int moduleIndex = -1;
        QTreeWidgetItem* top = sectionItem;
        while (top->parent()) top = top->parent();
        if (top->data(0, ModuleRole).isValid())
            moduleIndex = top->data(0, ModuleRole).toInt();
        if (moduleIndex >= 0 && moduleIndex < modules_.size())
            module = &modules_.at(moduleIndex);
    }

    if (module && !section.isEmpty()) {
        const QString lower = section.toLower();
        const bool blockMeshDocument = documentIndex >= 0
            && documentIndex < documents_.size()
            && documents_.at(documentIndex).fileName().compare(
                   QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0;
        const QVariant entryData = sectionItem->data(0, EntryRole);
        const bool isConcreteEntry = entryData.isValid();

        if (lower == QStringLiteral("vertices")) {
            if (isConcreteEntry && blockMeshDocument) {
                QAction* remove = menu.addAction(QStringLiteral("删除此顶点"));
                connect(remove, &QAction::triggered, this, [this, documentIndex, entryData] {
                    emit configEntryRemoveRequested(documentIndex, entryData.toInt());
                });
            } else {
                QAction* add = menu.addAction(QStringLiteral("新增顶点"));
                connect(add, &QAction::triggered, this, [this, documentIndex, section] {
                if (documentIndex >= 0 && documentIndex < documents_.size()
                    && documents_.at(documentIndex).fileName().compare(
                           QStringLiteral("blockMeshDict"), Qt::CaseInsensitive) == 0) {
                    emit configAssignmentsRequested(documentIndex, section,
                                                    {QStringLiteral("(0 0 0)")}, false);
                    return;
                }
                bool ok = false;
                const QString value = QInputDialog::getText(
                    this, QStringLiteral("新增顶点"), QStringLiteral("坐标 [x, y, z]"),
                    QLineEdit::Normal, QStringLiteral("[0, 0, 0]"), &ok);
                if (!ok || value.trimmed().isEmpty()) return;
                emit configAssignmentsRequested(documentIndex, section,
                    {QStringLiteral("v = %1").arg(value.trimmed())}, false);
                });
            }
            menu.addSeparator();
        }
        if (lower == QStringLiteral("blocks") && blockMeshDocument) {
            if (isConcreteEntry) {
                QAction* remove = menu.addAction(QStringLiteral("删除此 block"));
                connect(remove, &QAction::triggered, this, [this, documentIndex, entryData] {
                    emit configEntryRemoveRequested(documentIndex, entryData.toInt());
                });
            } else {
                QAction* add = menu.addAction(QStringLiteral("新增 block"));
                connect(add, &QAction::triggered, this, [this, documentIndex] {
                    emit configAssignmentsRequested(
                        documentIndex, QStringLiteral("blocks"),
                        {QStringLiteral("hex (0 0 0 0 0 0 0 0) (1 1 1) simpleGrading (1 1 1)")},
                        false);
                });
            }
            menu.addSeparator();
        }
        if (lower == QStringLiteral("boundary") && blockMeshDocument) {
            const QString patchName = sectionItem->data(0, PatchNameRole).toString();
            if (patchName.isEmpty()) {
                QAction* add = menu.addAction(QStringLiteral("新增 patch"));
                connect(add, &QAction::triggered, this, [this, documentIndex] {
                    QStringList names;
                    for (const ConfigEntry& entry : documents_.at(documentIndex).entries()) {
                        if (entry.section.startsWith(QStringLiteral("boundary/"),
                                                     Qt::CaseInsensitive)) {
                            names.push_back(entry.section.section(QLatin1Char('/'), -1));
                        }
                    }
                    int suffix = 0;
                    QString name;
                    do { name = QStringLiteral("patch%1").arg(suffix++); }
                    while (names.contains(name));
                    emit blockMeshPatchAddRequested(documentIndex, name,
                                                     QStringLiteral("patch"));
                });
            } else {
                QAction* addFace = menu.addAction(QStringLiteral("新增四边形面"));
                connect(addFace, &QAction::triggered, this, [this, documentIndex, patchName] {
                    emit blockMeshPatchFaceAddRequested(documentIndex, patchName,
                                                        QStringLiteral("(0 0 0 0)"));
                });
                QAction* remove = menu.addAction(QStringLiteral("删除此 patch"));
                connect(remove, &QAction::triggered, this, [this, documentIndex, patchName] {
                    emit blockMeshPatchRemoveRequested(documentIndex, patchName);
                });
            }
            menu.addSeparator();
        }
        if (lower == QStringLiteral("boundaryface") && blockMeshDocument
            && isConcreteEntry) {
            QAction* remove = menu.addAction(QStringLiteral("删除此面"));
            connect(remove, &QAction::triggered, this, [this, documentIndex, entryData] {
                emit configEntryRemoveRequested(documentIndex, entryData.toInt());
            });
            menu.addSeparator();
        }
    }

    if (module && module->kind == ProjectModuleKind::Multiphase) {
        QAction* phaseChange = menu.addAction(
            moduleIcon(ProjectModuleKind::PhaseChange),
            QStringLiteral("新建相变"));
        connect(phaseChange, &QAction::triggered, this, [this] {
            emit addModuleRequested(ProjectModuleKind::PhaseChange);
        });
        menu.addSeparator();
    }

    if (module && (module->kind == ProjectModuleKind::InitialCondition
                   || module->kind == ProjectModuleKind::BoundaryCondition)) {
        QStringList fields{QStringLiteral("p"), QStringLiteral("U")};
        if (module->densityBased) {
            fields << QStringLiteral("T") << QStringLiteral("rho")
                   << QStringLiteral("alpha");
        }
        for (const QString& field : fields) {
            QAction* addField = menu.addAction(
                QStringLiteral("新增字段 %1").arg(field));
            connect(addField, &QAction::triggered, this,
                    [this, kind = module->kind, field] {
                emit fieldAddRequested(kind, field);
            });
        }
        menu.addSeparator();
    }

    if (module && module->kind == ProjectModuleKind::Output) {
        const QString taskName = sectionItem
            ? sectionItem->data(0, OutputTaskRole).toString() : QString();
        if (taskName.isEmpty()) {
            QAction* create = menu.addAction(
                moduleIcon(ProjectModuleKind::Output), QStringLiteral("新建 job"));
            connect(create, &QAction::triggered,
                    this, &MainWindow::requestCreateOutputTask);
        } else {
            QAction* post = menu.addAction(
                moduleIcon(ProjectModuleKind::Output), QStringLiteral("查看后处理"));
            connect(post, &QAction::triggered, this, [this, taskName] {
                emit outputTaskSelected(taskName);
            });
        }
        menu.addSeparator();
    }

    addOptionalModuleActions(&menu);
    menu.exec(configTree_->viewport()->mapToGlobal(position));
}



void MainWindow::showOutputModule()
{
    focusModule(ProjectModuleKind::Output);
}


void MainWindow::handleTreeHover(QTreeWidgetItem* item)
{
    if (!item) { emit geometryHoverClear(); return; }
    QVariant blockVar = item->data(0, BlockIDRole);
    QVariant patchVar = item->data(0, PatchIDRole);
    QVariant faceVar  = item->data(0, FaceIDRole);
    if (blockVar.isValid())
        emit geometryHoverBlock(blockVar.toInt());
    else if (patchVar.isValid() && faceVar.isValid())
        emit geometryHoverFace(patchVar.toInt(), faceVar.toInt());
    else if (patchVar.isValid())
        emit geometryHoverPatch(patchVar.toInt());
    else
        emit geometryHoverClear();
}

} // namespace SF::GUI
