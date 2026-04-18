
# Chai-1

Chai-1 是一种先进的单链/复合物三维结构预测模型。在 ProtDesignTools 的封装下，它支持通过纯序列（FASTA）或直接提取 PDB 中的序列信息进行预测，并可以非常方便地进行特定链（如 DNA 或配体）的靶向替换。

## 支持的运行模式
- **序列预测**：从头基于字符串或 FASTA 文件预测结构。
- **复合物预测**：自动组装含有多条蛋白质链、DNA 链以及小分子配体（Ligand）的输入。
- **结构重构 (DNA/配体替换)**：输入现有 PDB 模板，保留蛋白质序列但替换目标 DNA 链。

## 命令行参数 (CLI Arguments)

### 核心输入参数 (三选一)
- `--pdb_path`: 结构文件路径，工具会自动提取其中的蛋白质链和 DNA 链。
- `--fasta_path`: 多链 FASTA 文件路径，如果提供则直接用作输入。
- `--sequence`: 单一蛋白质的序列字符串（例如 `"MASND..."`）。

### 复合体与配体参数
- `--target_dna_seq`: 当使用 `--pdb_path` 提取序列时，此参数用于将原本结构中的 DNA 链替换为你指定的目标 DNA 序列。
- `--ligand`: 字符串形式的配体 SMILES（如 `"CC(=O)OC1=CC=CC=C1C(=O)O"`）。工具会自动生成符合 Chai-1 解析规范的配体输入 Header。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认受全局配置控制，但强烈建议在 GPU 环境下运行)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/chai1/output`)。

---

## 命令行调用示例

### 1. 从纯序列预测单体结构
```bash
protdesign chai1 --sequence "MAQRTLEVW..."
```

### 2. 预测蛋白质与小分子配体的复合体
```bash
protdesign chai1 --sequence "MAQRTLEVW..." --ligand "CC(=O)OC1=CC=CC=C1C(=O)O"
```

### 3. 从现成 FASTA 文件预测多链复合体
```bash
protdesign chai1 --fasta_path ./data/complex.fasta
```

### 4. 从 PDB 提取序列并替换目标 DNA 链 (结构重构)
假设你有一个 DNA 结合蛋白的晶体结构，你想要预测它与一条全新 DNA 序列的结合模式：
```bash
protdesign chai1 --pdb_path ./data/dna_binding_protein.pdb --target_dna_seq "ATGCGTAC"
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回转换为 PDB 格式的预测结果以及解析好的置信度分数：

```python
from protdesigntools.tools.chai1.tool import Chai1

# 初始化工具
predictor = Chai1(exec_mode="slurm")

# 运行结构预测 (PDB提取 + DNA替换)
results = predictor(
    pdb_path="./data/dna_binding_protein.pdb",
    target_dna_seq="ATGCGTAC"
)

# 打印解析后的关键输出
print(f"Predicted Structure Path: {results['predicted_pdb']}")

# 置信度指标
print(f"pLDDT Score: {results['plddt']}")
print(f"iPTM Score: {results['iptm']}")
print(f"Combined Confidence Score: {results['combined_score']}")
```

### 输出结果说明
工具在执行完毕后，除了保存原始的 `.cif` 预测结果，还会自动完成以下几步操作（如果环境中安装了 `Biopython`）：
1. 自动定位分数最高（`model_idx_0`）的 `.cif` 文件。
2. 自动将其转换为标准的 `.pdb` 格式并返回路径（`predicted_pdb`）。
3. 扫描结构文件计算出平均 `pLDDT`。
4. 解析对应的 `.npz` 分数文件，提取 `iPTM`（针对复合物）。
5. 计算并返回 `combined_score` = `(pLDDT / 100.0 + iPTM) / 2.0`。
