/// @file SF_GeometryPanel.cpp
/// @brief GUI 几何编辑与参数面板实现。

#include "SF_GeometryPanel.h"
#include <QEvent>
#include <QHBoxLayout>

namespace SF::GUI {
namespace {

constexpr int kComboWidth = 42;
constexpr int kSeedSpinWidth = 56;
constexpr int kGradingSpinWidth = 56;

QDoubleSpinBox* makeSeedSpin(QWidget* p)
{
    auto* s = new QDoubleSpinBox(p);
    s->setRange(1, 99999); s->setDecimals(0);
    s->setFixedWidth(kSeedSpinWidth);
    s->setButtonSymbols(QAbstractSpinBox::NoButtons);
    s->setStyleSheet(QStringLiteral(
        "QDoubleSpinBox{background:#fffdf8;border:1px solid #aa9f8f;"
        "border-radius:4px;padding:1px 3px;color:#3d382f;font-size:11px;}"));
    return s;
}

QDoubleSpinBox* makeGradingSpin(QWidget* p)
{
    auto* s = new QDoubleSpinBox(p);
    s->setRange(0.01, 100.0); s->setDecimals(3);
    s->setSingleStep(0.1);
    s->setFixedWidth(kGradingSpinWidth);
    s->setButtonSymbols(QAbstractSpinBox::NoButtons);
    s->setStyleSheet(QStringLiteral(
        "QDoubleSpinBox{background:#fffdf8;border:1px solid #aa9f8f;"
        "border-radius:4px;padding:1px 3px;color:#3d382f;font-size:11px;}"));
    return s;
}

} // anonymous

// ============================================================
// CircleRadioGroup
// ============================================================

CircleRadioGroup::CircleRadioGroup(const QStringList& labels, QWidget* parent)
    : QWidget(parent)
{
    auto* row = new QHBoxLayout(this);
    row->setContentsMargins(0, 0, 0, 0);
    row->setSpacing(10);
    for (int i = 0; i < labels.size(); ++i) {
        auto* rb = new QRadioButton(labels.at(i), this);
        rb->setStyleSheet(QStringLiteral(
            "QRadioButton{spacing:4px;color:#4a443a;font-size:12px;}"
            "QRadioButton::indicator{width:14px;height:14px;border-radius:7px;"
            "border:1.6px solid #91846f;background:#fcf9f2;}"
            "QRadioButton::indicator:checked{background:qradialgradient("
            "cx:0.5,cy:0.5,radius:0.35,stop:0 #3d382f,stop:0.55 #3d382f,"
            "stop:0.65 #fcf9f2,stop:1 #fcf9f2);border:1.6px solid #3d382f;}"));
        connect(rb, &QRadioButton::toggled, this, [this,i](bool on){
            if (on) {
                for (int j = 0; j < buttons_.size(); ++j) buttons_[j]->setChecked(j==i);
                emit currentIndexChanged(i);
            }
        });
        buttons_.push_back(rb);
        row->addWidget(rb);
        if (i == 0) rb->setChecked(true);
    }
    row->addStretch();
}

int CircleRadioGroup::currentIndex() const
{
    for (int i = 0; i < buttons_.size(); ++i)
        if (buttons_.at(i)->isChecked()) return i;
    return 0;
}

void CircleRadioGroup::setCurrentIndex(int index)
{
    if (index >= 0 && index < buttons_.size())
        buttons_.at(index)->setChecked(true);
}

// ============================================================
// VerticesRow
// ============================================================

VerticesRow::VerticesRow(int pointCount, QWidget* parent) : QWidget(parent)
{
    auto* row = new QHBoxLayout(this);
    row->setContentsMargins(0,0,0,0); row->setSpacing(4);
    for (int i = 0; i < 4; ++i) {
        auto* cb = new QComboBox(this);
        cb->setFixedWidth(kComboWidth);
        cb->setStyleSheet(QStringLiteral(
            "QComboBox{background:#fffdf8;border:1px solid #aa9f8f;"
            "border-radius:4px;padding:1px 4px;color:#3d382f;font-size:11px;}"));
        cb->blockSignals(true);
        cb->clear();
        for (int j = 0; j < pointCount; ++j) cb->addItem(QString::number(j), j);
        cb->blockSignals(false);
        connect(cb, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &VerticesRow::valuesChanged);
        cb->installEventFilter(this);
        cb_[i] = cb;
        row->addWidget(cb);
    }
}

bool VerticesRow::eventFilter(QObject* watched, QEvent* event)
{
    for (QComboBox* combo : cb_) {
        if (watched != combo) continue;
        if (event->type() == QEvent::MouseButtonPress
            || event->type() == QEvent::FocusIn) {
            emit interactionStarted();
        }
        break;
    }
    return QWidget::eventFilter(watched, event);
}

QVector<int> VerticesRow::values() const
{
    QVector<int> v;
    for (int i = 0; i < 4; ++i) v.push_back(cb_[i]->currentData().toInt());
    return v;
}

void VerticesRow::setValues(int v0, int v1, int v2, int v3)
{
    const int vals[4] = {v0,v1,v2,v3};
    for (int i = 0; i < 4; ++i) {
        int idx = cb_[i]->findData(vals[i]);
        if (idx >= 0) cb_[i]->setCurrentIndex(idx);
    }
}

void VerticesRow::refreshPointRange(int pointCount)
{
    for (int i = 0; i < 4; ++i) {
        cb_[i]->blockSignals(true);
        int old = cb_[i]->currentIndex();
        cb_[i]->clear();
        for (int j = 0; j < pointCount; ++j) cb_[i]->addItem(QString::number(j), j);
        if (old >= 0 && old < pointCount) cb_[i]->setCurrentIndex(old);
        cb_[i]->blockSignals(false);
    }
}

// ============================================================
// TripleSpinRow
// ============================================================

TripleSpinRow::TripleSpinRow(const QString& label, bool integer, QWidget* parent)
    : QWidget(parent)
{
    auto* row = new QHBoxLayout(this);
    row->setContentsMargins(0,0,0,0); row->setSpacing(4);
    auto* lbl = new QLabel(label, this);
    lbl->setStyleSheet(QStringLiteral("color:#817767;font-size:11px;"));
    lbl->setFixedWidth(52);
    row->addWidget(lbl);
    for (int i = 0; i < 3; ++i) {
        auto* spin = integer ? makeSeedSpin(this) : makeGradingSpin(this);
        connect(spin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, &TripleSpinRow::valuesChanged);
        spins_[i] = spin;
        row->addWidget(spin);
    }
    row->addStretch();
}

QVector<double> TripleSpinRow::values() const
    { return {spins_[0]->value(), spins_[1]->value(), spins_[2]->value()}; }

void TripleSpinRow::setValues(double a, double b, double c)
    { spins_[0]->setValue(a); spins_[1]->setValue(b); spins_[2]->setValue(c); }

bool TripleSpinRow::isInteger() const { return spins_[0]->decimals() == 0; }

} // namespace SF::GUI
