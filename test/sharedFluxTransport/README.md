# Shared Flux Transport — Phase B

所有新测试、副本、日志、执行证据和运行临时目录都在本目录。不要把这些测试迁到 `/tmp`，不要覆盖正式 `test/Sod` / `test/IBM`。CFD endTime 不超过 1 s。

- `test_sharedFlux.cpp`：真实 NS face borrow；非均匀 rho 常浓度；开放边界守恒预算；Euler/RK4 时间和 Upwind1 空间收敛；1D/2D/3D 与变系数 Stage D。
- `test_stageSST.cpp`：原 k/omega 数组别名、保守 stage 输运、模型系数、真实 mu_t 黏性反馈、源/汇/扩散/对流隔离及能力负例。stage 路径中旧 correct() 调用即失败。
- `fluid_fixture.h`：沿用已有 scalar/flow harness 的配置工具，没有新的 runner。
- `run_short_tests.py`：读取当前 CTest inventory，迁移 native tests 的工作目录到本目录；运行 serial 子集。MPI 子集需要可启动 PRTE/socket 的 host 环境。
- `run_compatibility.py <新目录名>`：在新目录生成五对短时冻结副本并执行 before/after；已有目录拒绝覆盖。需要先保存当前任务冻结的 `evidence/sonicSolver-before`。
- `run_production.py`：准备 MPI Sod / native SST 副本和完整命令清单，**只准备，不自动执行**；已有目录拒绝覆盖。执行命令使用 manifest 中的绝对路径。
- `analyze_evidence.py`：最终 before/after VTS、dt 和带物理体积权重统计；排除 AppleDouble，MPI replica 不平均。

```sh
DEVELOPER_DIR=/Library/Developer/CommandLineTools cmake --build build --parallel 4
python3 tools/check_architecture.py
python3 test/sharedFluxTransport/run_short_tests.py
TMPDIR="$PWD/test/sharedFluxTransport/evidence" ctest --test-dir test/sharedFluxTransport/ctest-run --output-on-failure --parallel 1 -R 'distributed|haloFailureContract|peskinPartialRows|sparseCanonicalCopy'
python3 test/sharedFluxTransport/analyze_evidence.py
```

已有已完成证据请保留。`evidence/short-test-manifest.json` 明确列出未运行的原脚本，不能把子集写成全量 CTest 通过。generated case/result/log/binary 不提交；版本控制仅保留 C++、脚本和此说明。

精确网格、初值、边界、单位、recipe、最终时间/步数、误差与适用范围见 [中文完整报告](../../docs/shared-conservative-flux-transport-migration.md)。新 SST 使用显式 `models/turbulence.yaml: transportMethod: stageConservative`；旧 once-step/近壁 Ghost 结果不被新数学覆盖。新 scalar/SST MPI 和新 stage SST IBM/pressure 仍 Unsupported。
