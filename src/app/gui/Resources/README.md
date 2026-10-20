# GUI Resources

普通 UI 图标放在本目录下的 `icons/`：

- 树菜单图标，例如几何、网格、求解器、边界条件。
- 顶部模块下拉菜单图标。
- 工具按钮图标，例如运行、停止、导入、保存。
- 面板状态图标，例如启用、警告、错误。

启动动画帧不放在这里；它们属于 `../Splash/assets/`。

推荐格式是 SVG，建议画布使用 `24x24`，保持视觉主体落在 `20x20`
左右。PNG 只在需要位图质感时使用，常用尺寸是 `24x24` 和
`48x48(@2x)`。

当前 C++ 使用的是 `SF_resources.qrc` 里的 alias，例如：

- `:/sonic/icons/geometry.png`
- `:/sonic/icons/mesh.png`
- `:/sonic/icons/set-directory.png`

如果原始文件名以后要改，只要同步更新 qrc 的 `<file alias=...>`，
C++ 里的路径就不用跟着改。

macOS `.app` 在 Finder 和 Dock 里显示的图标使用
`icons/sonicGui.icns`，由 `icons/logo.png` 生成，并在
`Application/CMakeLists.txt` 里通过 `MACOSX_BUNDLE_ICON_FILE` 打包。
如果你换了 logo，重新生成这个 `.icns` 后再构建即可。
