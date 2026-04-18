
# ESM-2

ESM-2 是 Facebook (Meta) 开源的大规模蛋白质语言模型。在 ProtDesignTools 的封装下，它主要用于评估序列突变前后的 Fitness（打分），或者提取序列的高维表示（Representation）。

## 支持的运行模式
- **序列打分**：基于语言模型的似然性 (Likelihood) 评估序列突变后的适应度（Fitness）分数。

## 命令行参数 (CLI Arguments)

### 核心输入参数 (二选一)
- `--sequence`: 输入一条纯蛋白质序列字符串。
- `--fasta_path`: 输入包含单条或多条序列的 FASTA 文件。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/esm2/output`)。

---

## 命令行调用示例

### 1. 单序列打分
```bash
protdesign esm2 --sequence "MAQRTLEVW..."
```

### 2. 从 FASTA 批量打分
```bash
protdesign esm2 --fasta_path ./data/mutants.fasta
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回模型的输出分数。

```python
from protdesigntools.tools.esm2.tool import ESM2

# 初始化工具
scorer = ESM2(exec_mode="slurm")

# 运行序列打分任务
results = scorer(
    sequence="MAQRTLEVW..."
)

# 查看生成的打分结果
if results["status"] == "success":
    print(f"Sequence Score (Fitness): {results['sequence_score']}")
else:
    print("Scoring failed.")
```
