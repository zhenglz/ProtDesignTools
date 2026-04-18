
# OpenMM Simulation

OpenMM 是一个高性能的分子动力学（Molecular Dynamics, MD）模拟工具包。在 ProtDesignTools 的封装下，它主要用于对预测或设计的蛋白质结构进行物理构象驰豫（Relaxation），或进行较短时间的分子动力学模拟，以评估结构的动态稳定性和修正微观碰撞。

## 支持的运行模式
- **结构模拟**：输入 PDB 结构，通过指定的温度和步数进行隐式或显式溶剂环境下的构象优化和模拟。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--pdb_path`: 需要进行模拟或能量最小化的蛋白质三维结构 PDB 文件（必填）。

### 物理模拟参数
- `--temperature`: 模拟所用的绝对温度，单位为开尔文（Kelvin） (默认: `300` K)。
- `--steps`: MD 积分器（Integrator）运行的总时间步数 (默认: `10000` 步)。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定，建议由于计算量较大而使用 GPU 节点或集群)。
- `--output_dir`: 指定结果输出的文件夹路径 (默认: `work_dir/openmm/output`)。结果中将包含产生的轨迹文件（Trajectory, 默认 `.dcd`）和物理日志文件。

---

## 命令行调用示例

### 1. 基础结构模拟与能量最小化
在 300K 温度下对设计的突变体进行 50,000 步的动力学模拟，以检查结构是否会迅速散开：
```bash
protdesign openmm --pdb_path ./data/mutant_design.pdb --temperature 300 --steps 50000
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回模拟产生的轨迹和日志文件的绝对路径。

```python
from protdesigntools.tools.openmm.tool import OpenMMSimulation

# 初始化工具
md_runner = OpenMMSimulation(exec_mode="slurm")

# 运行动力学模拟任务
results = md_runner(
    pdb_path="./data/mutant_design.pdb",
    temperature=310.0,
    steps=20000
)

# 查看模拟输出文件
if results["status"] == "success":
    print(f"Trajectory file saved at: {results['trajectory']}")
    print(f"Simulation Log saved at: {results['log_file']}")
else:
    print("MD Simulation failed.")
```
