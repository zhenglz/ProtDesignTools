
# TMalign

TMalign 是一款经典的蛋白质三维结构比对工具。它可以基于序列无关（sequence-independent）的结构叠加算法，比较两个蛋白质结构的相似度，并输出 TM-score（一种衡量结构相似性的指标，不受局部较大偏差的影响）。

## 支持的运行模式
- **结构对齐**：比对靶标结构（Target）与参考结构（Reference），并输出对应的 RMSD、TM-score 和残基映射关系。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--target_pdb`: 待比对的目标 PDB 结构文件路径（必填）。
- `--reference_pdb`: 作为基准的参考 PDB 结构文件路径（必填）。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/tmalign/output`)。它会在此目录下生成对齐后的重叠坐标文件等详细数据。

---

## 命令行调用示例

### 1. 基础结构对齐
比对野生型（WT）与你新设计的突变体（Mutant）结构差异：
```bash
protdesign tmalign --target_pdb ./data/mutant.pdb --reference_pdb ./data/wildtype.pdb
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它能够返回一个结构化的结果字典，包含不同标准化基准下的 TM-score 和 RMSD。

```python
from protdesigntools.tools.tmalign.tool import TMalign

# 初始化工具
aligner = TMalign(exec_mode="local")

# 运行比对任务
results = aligner(
    target_pdb="./data/mutant.pdb",
    reference_pdb="./data/wildtype.pdb"
)

# 查看核心对齐指标
if results["status"] == "success":
    print(f"RMSD: {results['rmsd']} Å")
    print(f"TM-score (Normalized by Target): {results['tm_score_target']}")
    print(f"TM-score (Normalized by Reference): {results['tm_score_ref']}")
    
    # 获取序列级别的对应关系（残基映射）
    print("Alignment:")
    print(f"Ref: {results['alignment']['ref_seq']}")
    print(f"Tgt: {results['alignment']['target_seq']}")
else:
    print("Alignment failed.")
```
