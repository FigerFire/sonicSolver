# Icons

把你自己绘制的普通 GUI 图标放在这里。

命名建议使用小写短横线，例如：

- `module-geometry.svg`
- `module-mesh.svg`
- `module-solver.svg`
- `action-run.svg`
- `action-stop.svg`
- `status-warning.svg`

后续接入 Qt Resource 时，这些文件会通过类似
`:/sonic/icons/module-geometry.svg` 的路径在 C++ 里加载。

`sonicGui.icns` 是 macOS app bundle 图标，不走 Qt Resource。
它会被 CMake 复制到 `sonicGui.app/Contents/Resources/`。
