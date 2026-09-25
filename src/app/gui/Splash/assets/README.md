# Splash Assets

开屏动画帧放在 `frames/`，文件名必须按顺序编号：

- `frame-00.png`
- `frame-01.png`
- `frame-02.png`
- `frame-03.png`
- `frame-04.png`
- `frame-05.png`

`SF_splash.cpp` 只按顺序播放这些 alias，不关心图片来自哪个算例。
如果分发源码后想替换开屏，只需要覆盖这些 PNG，并保持文件名不变。

普通按钮、菜单和树节点图标不要放在这里；它们属于
`../../Resources/icons/`。
