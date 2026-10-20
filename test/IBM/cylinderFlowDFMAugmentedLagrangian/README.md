## 当前执行状态：Unsupported

本算例保留 cylinderFlowDFMAugmentedLagrangian 对应的原生 KKT 数学声明、无滑移约束和配置，没有改成其他 IBM 算法。

当前 provider 尚未实现它需要的 pressure predictor / block schedule；`sonicSolver check` 验证输入结构，可通过，但 `explain` 显示 Unsupported，实际 `run` 会在时间推进前 fail fast，未产生 CFD 解。

这不是缺一个 YAML 参数即可修复的问题。本轮只整理算例输入，没有补写 production KKT 数值内核。请先运行 `../cylinderFlowFictitiousDomain` 等已支持的方法；它们属于不同数学算法，不作为本方法的回归替代。
