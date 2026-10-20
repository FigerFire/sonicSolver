#pragma once

/// @file SF_GeometryPanel.h
/// @brief 配置树几何节点用的小组件：圆环单选、顶点下拉、三格数值。

#include <QComboBox>
#include <QDoubleSpinBox>
#include <QEvent>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QRadioButton>
#include <QString>
#include <QVector>
#include <QWidget>

namespace SF::GUI {

/// @brief 圆环风格的互斥单选组。选中时中心填充黑色圆点。
class CircleRadioGroup : public QWidget {
    Q_OBJECT
public:
    explicit CircleRadioGroup(const QStringList& labels, QWidget* parent = nullptr);
    int currentIndex() const;
    void setCurrentIndex(int index);
signals:
    void currentIndexChanged(int index);
private:
    QVector<QRadioButton*> buttons_;
};

/// @brief 一行四个顶点索引下拉框（block / patch face 共用）。
class VerticesRow : public QWidget {
    Q_OBJECT
public:
    explicit VerticesRow(int pointCount, QWidget* parent = nullptr);
    QVector<int> values() const;
    void setValues(int v0, int v1, int v2, int v3);
    void refreshPointRange(int pointCount);
    QComboBox* combo0() const { return cb_[0]; }
signals:
    void valuesChanged();
    /// @brief 用户开始操作某个顶点下拉框，用于让所属树节点成为当前选择。
    void interactionStarted();
private:
    bool eventFilter(QObject* watched, QEvent* event) override;
    QComboBox* cb_[4];
};

/// @brief 带标签的三个数值输入行（seed nx ny nz 或 grading gx gy gz）。
class TripleSpinRow : public QWidget {
    Q_OBJECT
public:
    explicit TripleSpinRow(const QString& label, bool integer, QWidget* parent = nullptr);
    QVector<double> values() const;
    void setValues(double a, double b, double c);
    bool isInteger() const;
signals:
    void valuesChanged();
private:
    QDoubleSpinBox* spins_[3];
};

} // namespace SF::GUI
