/// @file SF_mainWindowTree.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_mainWindow.h"

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

void MainWindow::rebuildTree() {
    rebuildingTree_=true;QSignalBlocker blocker(configTree_);configTree_->clear();
    for(int mi=0;mi<modules_.size();++mi) {
        const auto& module=modules_[mi];auto* root=new QTreeWidgetItem(configTree_);
        root->setText(0,module.title);root->setIcon(0,moduleIcon(module.kind));root->setData(0,ModuleRole,mi);
        if(module.kind==ProjectModuleKind::Output) {
            for(const auto& task:outputTasks_){auto* item=new QTreeWidgetItem(root);item->setText(0,task.name);item->setText(1,displayTime(task)+" / "+displayProgress(task.progress));item->setData(0,OutputTaskRole,task.name);item->setToolTip(0,task.pvdPath);}
        } else for(const auto& ref:module.entries) {
            const auto& doc=documents_.at(ref.documentIndex);const auto& e=doc.entries().at(ref.entryIndex);
            auto* item=new QTreeWidgetItem(root);item->setText(0,e.key);item->setText(1,"Edit registered native YAML…");item->setToolTip(1,e.value);
            item->setData(0,DocumentRole,ref.documentIndex);item->setData(0,EntryRole,ref.entryIndex);
        }
        root->setExpanded(mi==currentModuleIndex_);
    }
    rebuildingTree_=false;updateEnabledMenu();
}

QString MainWindow::entryLabel(const ConfigEntry& entry) const {return entry.key;}



void MainWindow::editTreeItem(QTreeWidgetItem* item) {
    if(!item || handleOutputTreeClick(item) || !item->data(0,DocumentRole).isValid())return;
    const int di=item->data(0,DocumentRole).toInt(),ei=item->data(0,EntryRole).toInt();
    bool ok=false;const auto& entry=documents_.at(di).entries().at(ei);
    const auto text=QInputDialog::getMultiLineText(this,"Native YAML",documents_.at(di).path(),entry.value,&ok);
    if(ok && text!=entry.value)emit configValueEdited(di,ei,text);
}

void MainWindow::commitTreeItemEdit(QTreeWidgetItem*,int) {}

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

QStringList MainWindow::enumChoices(const ConfigEntry&) const {return {};}

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
    for(const auto& module:modules_)if(module.optional && module.enabled && module.documentIndex>=0) {
        auto* action=enabledMenu_->addAction(module.title);
        connect(action,&QAction::triggered,this,[this,index=module.documentIndex]{selectEntry(index,0);});
    }
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

void MainWindow::showModuleContextMenu(const QPoint& position) {
    QMenu menu;auto* item=configTree_->itemAt(position);
    if(item && item->data(0,DocumentRole).isValid()) {
        auto* edit=menu.addAction("Edit registered native YAML");connect(edit,&QAction::triggered,this,[this,item]{editTreeItem(item);});
    }
    for(const auto& entry:QVector<QPair<QString,ProjectModuleKind>>{{"Register turbulence (inactive)",ProjectModuleKind::Turbulence},{"Register IBM (inactive)",ProjectModuleKind::IBM},{"Register MRF (inactive)",ProjectModuleKind::MRF}}) {
        auto* add=menu.addAction(entry.first);connect(add,&QAction::triggered,this,[this,kind=entry.second]{emit addModuleRequested(kind);});
    }
    menu.exec(configTree_->viewport()->mapToGlobal(position));
}



void MainWindow::showOutputModule()
{
    focusModule(ProjectModuleKind::Output);
}


void MainWindow::handleTreeHover(QTreeWidgetItem*) {emit geometryHoverClear();}

} // namespace SF::GUI
