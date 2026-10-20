/// @file SF_splash.h
/// @brief GUI 启动画面资源与动画控制。

#pragma once

#include <QStringList>
#include <QWidget>

class QLabel;
class QTimer;

namespace SF::GUI {

/// @brief 应用启动阶段显示的压力切片动画窗口。
class SplashScreen final : public QWidget {
    Q_OBJECT

public:
    explicit SplashScreen(QWidget* parent = nullptr);

    /// @brief 将启动窗口放到当前屏幕中央并显示。
    void showCentered();

    /// @brief 更新当前启动阶段文字。
    void setStage(const QString& stage);

    /// @brief 淡出启动窗口，并激活主窗口。
    void finish(QWidget* mainWindow);

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    void advanceFrame();
    void updateFrame();

    QStringList frames_;
    int frameIndex_ = 0;
    QLabel* imageLabel_ = nullptr;
    QLabel* captionLabel_ = nullptr;
    QLabel* stageLabel_ = nullptr;
    QTimer* frameTimer_ = nullptr;
};

} // namespace SF::GUI
