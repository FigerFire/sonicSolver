/// @file SF_splash.cpp
/// @brief GUI 启动画面资源与动画控制。

#include "SF_splash.h"

#include <QApplication>
#include <QCursor>
#include <QGraphicsDropShadowEffect>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainter>
#include <QPixmap>
#include <QPropertyAnimation>
#include <QScreen>
#include <QTimer>
#include <QVBoxLayout>

static void initializeSplashResources()
{
    Q_INIT_RESOURCE(SF_splash);
    Q_INIT_RESOURCE(SF_resources);
}

static QStringList splashFramePaths()
{
    QStringList frames;
    for (int i = 0; i < 6; ++i) {
        frames.push_back(QStringLiteral(":/sonic/splash/frame-%1.png")
                             .arg(i, 2, 10, QLatin1Char('0')));
    }
    return frames;
}

namespace SF::GUI {

SplashScreen::SplashScreen(QWidget* parent)
    : QWidget(parent)
    , frames_(splashFramePaths())
{
    initializeSplashResources();
    setWindowFlags(Qt::SplashScreen | Qt::FramelessWindowHint
                   | Qt::WindowStaysOnTopHint);
    setAttribute(Qt::WA_TranslucentBackground);
    setFixedSize(780, 420);

    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(34, 30, 34, 30);
    root->setSpacing(12);

    auto* titleRow = new QWidget(this);
    auto* titleLayout = new QHBoxLayout(titleRow);
    titleLayout->setContentsMargins(0, 0, 0, 0);
    titleLayout->setSpacing(14);

    auto* logoLabel = new QLabel(titleRow);
    logoLabel->setFixedSize(58, 58);
    logoLabel->setAlignment(Qt::AlignCenter);
    const QPixmap logo(QStringLiteral(":/sonic/icons/logo.png"));
    logoLabel->setPixmap(logo.scaled(
        QSize(50, 50), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    logoLabel->setStyleSheet(QStringLiteral(
        "background:rgba(255,252,246,185);"
        "border:1px solid rgba(145,132,111,80);"
        "border-radius:15px;"));
    titleLayout->addWidget(logoLabel);

    auto* titleText = new QWidget(titleRow);
    auto* titleTextLayout = new QVBoxLayout(titleText);
    titleTextLayout->setContentsMargins(0, 2, 0, 0);
    titleTextLayout->setSpacing(2);

    auto* title = new QLabel(QStringLiteral("sonicSolver"), titleText);
    QFont titleFont = title->font();
    titleFont.setPointSizeF(25);
    titleFont.setBold(true);
    title->setFont(titleFont);
    title->setStyleSheet(QStringLiteral("color:#3d382f;"));
    titleTextLayout->addWidget(title);

    auto* subtitle = new QLabel(
        QStringLiteral("Structured CFD Workbench"), titleText);
    subtitle->setStyleSheet(
        QStringLiteral("color:#817767; letter-spacing:1px;"));
    titleTextLayout->addWidget(subtitle);
    titleTextLayout->addStretch();
    titleLayout->addWidget(titleText, 1);
    root->addWidget(titleRow);

    imageLabel_ = new QLabel(this);
    imageLabel_->setFixedHeight(190);
    imageLabel_->setAlignment(Qt::AlignCenter);
    imageLabel_->setStyleSheet(QStringLiteral(
        "background:rgba(255,252,246,205);"
        "border:1px solid rgba(145,132,111,95);"
        "border-radius:14px;"));
    root->addWidget(imageLabel_);

    captionLabel_ = new QLabel(this);
    captionLabel_->setStyleSheet(
        QStringLiteral("color:#5f574b; font-weight:600;"));
    root->addWidget(captionLabel_);

    stageLabel_ = new QLabel(QStringLiteral("正在初始化界面..."), this);
    stageLabel_->setStyleSheet(QStringLiteral("color:#817767;"));
    root->addWidget(stageLabel_);
    root->addStretch();

    auto* shadow = new QGraphicsDropShadowEffect(this);
    shadow->setBlurRadius(38);
    shadow->setOffset(0, 12);
    shadow->setColor(QColor(72, 63, 50, 90));
    imageLabel_->setGraphicsEffect(shadow);

    frameTimer_ = new QTimer(this);
    frameTimer_->setInterval(330);
    connect(frameTimer_, &QTimer::timeout,
            this, &SplashScreen::advanceFrame);
    updateFrame();
}

void SplashScreen::showCentered()
{
    QScreen* screen = QGuiApplication::screenAt(QCursor::pos());
    if (!screen) screen = QGuiApplication::primaryScreen();
    if (screen) {
        const QRect available = screen->availableGeometry();
        move(available.center() - rect().center());
    }
    setWindowOpacity(1.0);
    show();
    raise();
    frameTimer_->start();
    qApp->processEvents();
}

void SplashScreen::setStage(const QString& stage)
{
    stageLabel_->setText(stage);
    advanceFrame();
    qApp->processEvents();
}

void SplashScreen::finish(QWidget* mainWindow)
{
    frameTimer_->stop();
    if (mainWindow) {
        mainWindow->show();
        mainWindow->raise();
        mainWindow->activateWindow();
    }
    auto* fade = new QPropertyAnimation(this, "windowOpacity", this);
    fade->setDuration(240);
    fade->setStartValue(1.0);
    fade->setEndValue(0.0);
    connect(fade, &QPropertyAnimation::finished, this, &QWidget::hide);
    fade->start(QAbstractAnimation::DeleteWhenStopped);
}

void SplashScreen::paintEvent(QPaintEvent*)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    const QRectF card = rect().adjusted(12, 10, -12, -18);
    painter.setPen(QPen(QColor(145, 132, 111, 95), 1.0));
    painter.setBrush(QColor(244, 238, 228, 244));
    painter.drawRoundedRect(card, 26, 26);
}

void SplashScreen::advanceFrame()
{
    if (frames_.isEmpty()) return;
    frameIndex_ = (frameIndex_ + 1) % frames_.size();
    updateFrame();
}

void SplashScreen::updateFrame()
{
    if (frames_.isEmpty()) return;
    const QPixmap frame(frames_.at(frameIndex_));
    imageLabel_->setPixmap(frame.scaled(
        imageLabel_->size() - QSize(20, 20),
        Qt::KeepAspectRatio, Qt::SmoothTransformation));
    captionLabel_->setText(QStringLiteral("CFD pressure-field preview"));
}

} // namespace SF::GUI
