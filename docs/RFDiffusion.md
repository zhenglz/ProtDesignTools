
# RFDiffusion (RFD3)

RFDiffusion (基于 RoseTTAFold Diffusion 的第三代工具) 是一款强大的蛋白质结构生成和骨架设计模型。在 ProtDesignTools 的封装下，它主要用于结构修复、从头设计（De novo design）和片段组装（Scaffolding）。

## 支持的运行模式
- **基于模板的设计**：输入已知 PDB 和 `contig` (结构约束字符串)，生成补全后的 `.cif.gz` 或 `.pdb` 结构。
- **从头生成**：在不提供结构模板的情况下，仅通过长度约束生成全新骨架。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--pdb_path`: 结构模板文件路径（必填，如果使用基于模板的设计）。
- `--contig`: 定义片段约束的字符串（例如 `"10-100"`，表示在模板中留出 10 到 100 个氨基酸长度的空隙用于模型重新设计骨架）。如果不提供，系统将尝试生成默认的无条件（Unconditional）骨架。
- `--length`: 限制生成骨架的总长度范围（例如 `"50-100"`）。
- `--fixed_atoms`: JSON 字符串，用于指定必须在空间坐标中固定的特定原子。
- `--num_designs`: 需要生成的设计数量 (默认: `1`)。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/rfdiffusion/output`)。

---

## 命令行调用示例

### 1. 基础骨架补全 (Scaffolding)
假设你的输入 PDB 在中间有缺失，你想让模型自动生成一段长度在 10-50 之间的连结序列（Linker）：
```bash
protdesign rfdiffusion --pdb_path ./data/scaffold.pdb --contig "A1-10/10-50/A60-100" --num_designs 5
```

### 2. 从头骨架生成 (Unconditional Design)
生成一段长度在 50-100 氨基酸之间的全新结构骨架：
```bash
protdesign rfdiffusion --pdb_path ./data/dummy.pdb --length "50-100" --contig "50-100"
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回生成的所有压缩结构文件的路径列表：

```python
from protdesigntools.tools.rfdiffusion.tool import RFDiffusion

# 初始化工具
rfd3 = RFDiffusion(exec_mode="slurm")

# 运行设计任务
results = rfd3(
    pdb_path="./data/scaffold.pdb",
    contig="10-100",
    num_designs=3
)

# 查看生成的结构文件路径
if results["status"] == "success":
    for idx, filepath in enumerate(results["generated_files"]):
        print(f"Design {idx + 1} saved at: {filepath}")
else:
    print("Design failed.")
```

### 输出结果说明
工具在执行时，会自动在临时目录生成用于推理的 `inference.json` 配置文件。执行完毕后，所有生成的 `.cif.gz` 或 `.pdb` 设计文件将被移动到你指定的 `output_dir` 中。`generated_files` 键中包含所有输出文件的绝对路径列表。
