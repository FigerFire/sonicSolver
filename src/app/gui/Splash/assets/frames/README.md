# Splash Frames

这里的 PNG 按文件名顺序播放。建议保持同一宽高比例，当前 splash 预览框
会自动按比例缩放图片。

替换规则：

1. 保持文件名 `frame-00.png` 到 `frame-05.png` 不变。
2. 保持透明或实底都可以，Qt 会缩放到开屏预览框。
3. 如果要增加或减少帧数，需要同步修改 `SF_splash.qrc` 和
   `splashFramePaths()`。
