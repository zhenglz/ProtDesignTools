
# DLKcat

DLKcat 是一个基于深度学习的酶催化活性（kcat）预测模型。它能够根据输入的酶氨基酸序列和底物的 SMILES 结构，预测酶对该底物的周转数（Turnover number）。

## 支持的运行模式
- **kcat 预测**：给定酶序列和底物化学结构，预测 $k_{cat}$ 值。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--sequence`: 酶的纯氨基酸序列字符串（必填）。
- `--smiles`: 底物的化学结构表示，采用标准的 SMILES 字符串格式（必填）。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/dlkcat/output`)。

---

## 命令行调用示例

### 1. 基础活性预测
假设你要预测某条序列针对阿司匹林（Aspirin, `CC(=O)OC1=CC=CC=C1C(=O)O`）的催化活性：
```bash
protdesign dlkcat --sequence "MAQRTLEVW..." --smiles "CC(=O)OC1=CC=CC=C1C(=O)O"
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回包含预测数值的字典：

```python
from protdesigntools.tools.dlkcat.tool import DLKcat

# 初始化工具
predictor = DLKcat(exec_mode="slurm")

# 运行活性预测任务
results = predictor(
    sequence="MAQRTLEVW...",
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O"
)

# 查看生成的预测结果
if results["status"] == "success":
    print(f"Predicted kcat: {results['kcat_prediction']} {results['unit']}")
else:
    print("Prediction failed.")
```
