
# Pythia

Pythia 是一种专用于评估蛋白质热稳定性变化的深度学习模型。它通过评估给定突变（Mutation）前后结构的差异，快速估算出自由能变化（$\Delta\Delta G$ 或 ddG），从而判断突变是否有利于蛋白质的折叠稳定性。

## 支持的运行模式
- **热稳定性评估**：给定参考 PDB 和一组点突变序列定义，输出计算得到的 ddG 值。通常，负数代表突变是稳定的（Stabilizing），正数代表不稳定。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--pdb_path`: 野生型（Wild Type）参考 PDB 文件路径（必填）。
- `--mutations`: 逗号分隔的点突变字符串，如 `"A12G, S30C"`（必填）。该格式包含原氨基酸、位置及突变后的氨基酸。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/pythia/output`)。

---

## 命令行调用示例

### 1. 基础单点突变打分
评估野生型蛋白在第 198 位由蛋氨酸（M）突变为苯丙氨酸（F）带来的热稳定性变化：
```bash
protdesign pythia --pdb_path ./data/wildtype.pdb --mutations "M198F"
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以返回具体打分的字典对象。

```python
from protdesigntools.tools.pythia.tool import Pythia

# 初始化工具
ddg_scorer = Pythia(exec_mode="slurm")

# 运行稳定性预测任务
results = ddg_scorer(
    pdb_path="./data/wildtype.pdb",
    mutations="M198F"
)

# 查看核心预测指标
if results["status"] == "success":
    print(f"Evaluated Mutation(s): {results['mutations']}")
    print(f"Predicted ddG: {results['ddg']} kcal/mol")
    
    if results['ddg'] < 0:
        print("Conclusion: The mutation is likely STABILIZING.")
    else:
        print("Conclusion: The mutation is likely DESTABILIZING.")
else:
    print("Pythia ddG prediction failed.")
```
