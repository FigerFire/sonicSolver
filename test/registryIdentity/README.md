# Registry identity regression

本目录验证 native IO 的母子身份和内部角色。正式规则见 AGENTS.md §34 与 src/ARCHITECTURE.md。

- `test/test_registryIdentity.cpp` 是自动回归：子标题、保留键/角色、未注册文件、目录/重复引用和 symlink。
- `migrate_cases.py` 是显式的一次性旧 fixture 转换工具，不是 production reader，不会由 solver 自动调用。调用前保留输入快照；遇到错误即停止，需检查已转换输入，不使用旧格式 fallback。
- `evidence/`、`*-before/`、`*-after/`、`ctest-run/`、`cases/` 都是本地验证数据，不上传 Git。before 保存本轮开始时的 executable/配置，after 只改变 IO 格式；全部求解输出位于对应 case/result/。
- 验证配置、命令、精确结果及 Unsupported 说明见 `docs/registry-identity-io-closure.md`。本轮不把历史 level-set 脚本的既有 capability 拒绝解释成 IO 或数值迁移成功。

运行编译后的 `./build/test_registryIdentity`；可传入一个 case 目录进行仅 CaseIO 的注册读取检查。数值回归必须使用当前 native 文档，旧冻结输入由旧 binary 读取，不能混用。
