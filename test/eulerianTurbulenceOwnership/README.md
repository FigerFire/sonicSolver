# Phase C：Eulerian turbulence ownership 验证

本目录只验证当前工作区的职责迁移，不替代历史 production baseline、长期湍流物理或新增 MPI 能力验收。正式算例不被改写。

所有 CFD startTime=0，Eulerian endTime=0.0002 s，保护案例 endTime=0.001 s（固定密度 PISO 为 0.0002 s）。嵌入式 scalar harness 的最长 endTime=1 s。每个 case 的结果写入自己的 `result/`。MPI 单测的 TMPDIR 与工作目录均在本目录内。

## 证据

- `evidence/sonicSolver-before`：迁移前当前工作区可执行文件。
- `baseline/`：10 个原执行文件的短时配置与结果。
- `instrumented/`：迁移前加入只读诊断的程序；先证明输出与 baseline 完全一致。
- `c1/`、`c2/`、`c3/`、`final/`：每阶段同输入运行；比较物理输出、矩阵、解与诊断，不放宽 tolerance。
- `protection/<name>-before|after/`：Sod、Ghost、stage SST、legacy SST/Ghost、固定密度 PISO、Peskin、FTS、fractional DLM、DFM 的隔离副本。
- `evidence/metrics.json`：每帧场 min/max、L1 mean、L2 RMS、Linf、Cartesian nodal volume integral，以及 old-time/矩阵/线性解/trace 统计。
- `evidence/case-settings.md`：从实际注册输入和实际日志生成的逐例设置。包括 WHAT、STATE、EOS、模型、边界、HOW/WHICH、网格、运行时间和路径，不能只看 case 名称推测。
- `evidence/*-manifest.json`、`executable-identities.json`：命令、退出码、输入/执行文件标识、严格对照结果。
- 中文结论：`docs/eulerian-turbulence-execution-ownership-migration.md`。

## 可重复验证

从仓库根目录执行；需要已有冻结的 before/instrumented 证据。`setup.py`、`freeze.py` 只用于迁移前建基线，已有冻结程序时会拒绝覆盖。不能在新代码上重建并冒充旧基线。

```sh
DEVELOPER_DIR=/Library/Developer/CommandLineTools cmake --build build --parallel 4
python3 tools/check_architecture.py
python3 test/eulerianTurbulenceOwnership/compare.py final
python3 test/eulerianTurbulenceOwnership/run_protection.py
python3 test/eulerianTurbulenceOwnership/analyze.py
python3 test/eulerianTurbulenceOwnership/settings.py
python3 test/eulerianTurbulenceOwnership/run_short_tests.py
TMPDIR="$PWD/test/eulerianTurbulenceOwnership/evidence" ctest \
  --test-dir test/eulerianTurbulenceOwnership/ctest-run --output-on-failure --parallel 1 \
  -R '^(distributed|haloFailureContract|peskinPartialRows|sparseCanonicalCopy)'
```

最后一项需 host MPI socket 权限。选中的 40 serial＋10 MPI/HYPRE C++ tests 不是根 build 的完整 67-test suite；其余脚本可能使用旧 IO、长时或临时目录，未在本轮执行。排除原因记录在 `short-test-manifest.json`。

## 单元与负例

`test_ownership.cpp` 的独立 phase harness 使用两个不同的 phaseMass、反向相选择和非均匀 scalar 数组。验证 physical alias、所有 ghost layers 的原边界 stencil、phase slot、缺失/别名 old-time、系数 publication/epoch、非法 mass/sink/系数、失败后 old-time 不变和重复 commit。`referenceBoundary.h` 是测试内的原 stencil 对照，不进入生产 library。

`test/test_eulerianTurbulence.cpp` 验证原生数学/STATE/HOW/WHICH contract、模型/用户等价、OpRegistry owner、三层循环及新增 solve 尾部顺序，并覆盖缺状态、错 backing、缺 pair、重复 occurrence/commit、未知/重复相、不支持的梯度/输运 order、MPI capability。

模型、物理数组和 closure 仍有真实职责；本轮只移动 solver workspace/BC/assembly/solve/commit，不用单流体显式 SST 覆盖 Eulerian 相质量隐式输运。
