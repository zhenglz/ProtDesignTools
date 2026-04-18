
# ProteinMPNN / LigandMPNN

ProteinMPNN 及其衍生模型（如 LigandMPNN）是一种基于图神经网络的序列设计和打分工具。它能够在给定蛋白质主链结构（或特定配体上下文）的情况下，生成具有高稳定性的氨基酸序列，或对现有序列和突变体进行打分评估。

## 支持的运行模式
- **设计模式 (design)**：基于 PDB 结构生成新的氨基酸序列。
- **打分模式 (scoring)**：对输入的结构和序列进行一致性打分。

## 命令行参数 (CLI Arguments)

### 通用参数
- `--mode`: 运行模式，可选 `design` 或 `scoring` (默认: `design`)。
- `--pdb_path`: 核心输入参数，必须提供一个 PDB 文件路径作为结构模板。
- `--exec_mode`: 运行方式，可选 `local` (本地执行) 或 `slurm` (提交到集群) (默认: 全局配置决定)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/proteinmpnn/output`)。

### 设计模式专属参数
- `--model_type`: 选择底层的模型版本，支持 `protein_mpnn`, `ligand_mpnn`, `soluble_mpnn`, `global_mpnn` (默认: `protein_mpnn`)。
- `--num_seqs`: 每个目标生成的序列数量 (默认: `1`)。
- `--sampling_temp`: 采样温度，数值越小越趋近于保守模型解，越大则增加序列多样性 (默认: `0.1`)。
- `--design_chains`: 指定需要设计的链，多条链用逗号分隔，如 `"A,B"`。未指定的链将作为固定的上下文存在。
- `--redesigned_residues`: (仅 LigandMPNN 支持) 显式指定需要重新设计的残基位置，如 `"A12 A13 B2"`。
- `--pack_side_chains`: (仅 LigandMPNN 支持) 是否打包侧链构象，`1` 为开启，`0` 为关闭。
- `--pack_with_ligand_context`: (仅 LigandMPNN 支持) 打包侧链时是否考虑配体环境，`1` 为开启。

#### 高级位置限制控制 (Classic ProteinMPNN 支持)
- `--fixed_positions`: JSON 字符串，锁定特定链的残基位置不被突变。如 `'{"A": [1, 2, 3]}'`。
- `--omit_AAs`: 全局忽略生成特定的氨基酸。如 `"CX"`（不生成半胱氨酸和未知残基）。
- `--omit_AA_per_pos`: JSON 字符串，在特定位置忽略特定氨基酸。如 `'{"A": {"1": "C", "2": "FWY"}}'`。
- `--bias_AA_per_pos`: JSON 字符串，为特定位置的特定氨基酸添加能量偏好值。如 `'{"A": {"1": {"A": 10.0}}}'`。
- `--tied_positions`: JSON 字符串，强制不同链的特定位置生成相同的氨基酸（适用于对称同源多聚体）。如 `'{"A": [1,2], "B": [1,2]}'`。

### 打分模式专属参数
- `--mutations`: 逗号分隔的突变字符串，如 `"A12G, S30C"`。如果提供，工具会自动将 PDB 提取出的序列应用突变后再进行打分。
- `--fasta_path`: 绕过 PDB 直接输入现成的 FASTA 文件进行打分。
- `--sequence`: 直接输入一条氨基酸字符串进行打分。

---

## 命令行调用示例

### 1. 基础设计 (生成 5 条新序列)
```bash
protdesign proteinmpnn --mode design --pdb_path ./data/example.pdb --num_seqs 5
```

### 2. 使用 LigandMPNN 进行侧链打包与设计
```bash
protdesign proteinmpnn --mode design --model_type ligand_mpnn --pdb_path ./data/complex.pdb --design_chains A --pack_side_chains 1
```

### 3. 高级限制设计 (固定 A 链的 1,2,3 位，且全局不生成半胱氨酸)
```bash
protdesign proteinmpnn --mode design --pdb_path ./data/example.pdb --fixed_positions '{"A": [1, 2, 3]}' --omit_AAs C
```

### 4. 突变体打分
```bash
protdesign proteinmpnn --mode scoring --pdb_path ./data/example.pdb --mutations "A12G, S30C"
```

---

## Python API 调用示例 (Import)

作为 Python 包集成到你的脚本中，能够直接返回包含解析结果的字典：

```python
from protdesigntools.tools.proteinmpnn.tool import ProteinMPNN

# 初始化工具 (可覆盖执行模式)
mpnn = ProteinMPNN(exec_mode="slurm")

# 1. 运行设计模式
design_results = mpnn(
    mode="design",
    model_type="ligand_mpnn",
    pdb_path="./data/example.pdb",
    design_chains="A",
    num_seqs=10,
    sampling_temp=0.1
)

# 查看解析出的设计序列及其分数
for seq_info in design_results["sequences"]:
    print(f"Sequence: {seq_info['seq']}")
    print(f"Score: {seq_info['score']}")
    print(f"Recovery: {seq_info.get('seq_recovery', 'N/A')}")

# 2. 运行突变打分模式
scoring_results = mpnn(
    mode="scoring",
    pdb_path="./data/example.pdb",
    mutations="M198F"
)

print(f"Mutant Score: {scoring_results['score']}")
```
