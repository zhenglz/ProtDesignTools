
# AutoDock Vina

AutoDock Vina 是一个著名的开源分子对接（Molecular Docking）程序。在 ProtDesignTools 的封装下，它主要用于预测小分子配体（Ligand）与受体蛋白质（Receptor）结合位点及亲和力的初步估算，能够输出结合能分数并生成候选姿态。

## 支持的运行模式
- **分子对接**：给定受体结构、配体结构以及包含结合位点坐标的三维搜索网格（Grid Box），执行刚性受体-柔性配体的构象搜索。

## 命令行参数 (CLI Arguments)

### 核心输入参数
- `--receptor`: 蛋白质受体文件路径（PDB 或 PDBQT 格式，必填）。
- `--ligand`: 小分子配体文件路径（SDF 或 PDBQT 格式，必填）。

### 搜索网格 (Grid Box) 参数
- `--center_x`: 对接搜索框中心点的 X 坐标 (默认: `0.0`)。
- `--center_y`: 对接搜索框中心点的 Y 坐标 (默认: `0.0`)。
- `--center_z`: 对接搜索框中心点的 Z 坐标 (默认: `0.0`)。
- `--size_x`: 搜索框在 X 轴上的尺寸 (默认: `20.0` Å)。
- `--size_y`: 搜索框在 Y 轴上的尺寸 (默认: `20.0` Å)。
- `--size_z`: 搜索框在 Z 轴上的尺寸 (默认: `20.0` Å)。

### 通用参数
- `--exec_mode`: 运行方式，可选 `local` 或 `slurm` (默认: 全局配置决定)。
- `--output_dir`: 指定对接结果输出的文件夹路径 (默认: `work_dir/vina/output`)。结果中将包含一系列结合姿态的 PDBQT 文件及相应的亲和力日志。

---

## 命令行调用示例

### 1. 基础配体对接
在一个 20x20x20 埃的网格内进行小分子与蛋白质活性中心的对接搜索：
```bash
protdesign vina --receptor ./data/protein.pdbqt --ligand ./data/aspirin.pdbqt --center_x 12.5 --center_y -4.2 --center_z 8.0 --size_x 20 --size_y 20 --size_z 20
```

---

## Python API 调用示例 (Import)

作为 Python 模块集成，它可以直接返回最佳结合能。

```python
from protdesigntools.tools.autodock_vina.tool import AutoDockVina

# 初始化工具
docker = AutoDockVina(exec_mode="local")

# 运行对接任务
results = docker(
    receptor="./data/protein.pdbqt",
    ligand="./data/aspirin.pdbqt",
    center_x=12.5,
    center_y=-4.2,
    center_z=8.0,
    size_x=20,
    size_y=20,
    size_z=20
)

# 查看对接核心指标
if results["status"] == "success":
    print(f"Top Binding Energy: {results['binding_energy']} kcal/mol")
    print(f"Docked Poses Saved At: {results['docked_file']}")
else:
    print("Docking failed.")
```
