
# ProtDesignTools

这是一个模块化、标准化的蛋白质设计工具集合，支持本地（Local）和集群（Slurm）调用范式。

## 核心设计理念

1.  **积木式调用**：通过定义标准 I/O（JSON）和调用接口，实现工具的灵活组合。
2.  **多模式支持**：同一工具可支持不同工作模式（如生成、打分、模拟）。
3.  **任务管理**：内置 `TaskManager`，支持 Slurm 状态追踪、最大任务数控制及全生命周期管理。
4.  **灵活调用**：支持 Python 脚本调用（CLI）和 `import` 调用（API）。

## 项目结构

- `core/`: 核心基类
  - `task_manager.py`: 任务管理（Slurm/Local）
  - `sequence.py`: 序列与突变处理
  - `structure.py`: 结构与多链管理（继承自 Sequence）
  - `base_tool.py`: 工具基类，定义标准 I/O 范式
- `tools/`: 具体工具实现（如 ProteinMPNN, AF3, Chai-1 等）
- `protocols/`: 流程化协议示例
- `data/`: 配置文件与示例数据

## 快速开始

### 1. CLI 调用方式

```bash
python main.py proteinmpnn --config tools/proteinmpnn/config.json
```

### 2. Python API 调用方式

```python
from tools.proteinmpnn.tool import ProteinMPNN

tool = ProteinMPNN()
results = tool.run({
    "mode": "design",
    "pdb_path": "example.pdb",
    "num_seqs": 5,
    "exec_mode": "slurm"
})
print(results)
```

## 扩展新工具

1. 在 `tools/` 下创建新文件夹。
2. 继承 `BaseTool` 类并实现 `run` 方法。
3. 定义该工具的标准 JSON 输入格式。
4. 在 `main.py` 中注册新工具。

## 任务管理说明

- **本地模式**：使用 `subprocess` 管理进程。
- **Slurm 模式**：自动生成 `.sh` 脚本，提交并追踪状态（PENDING -> RUNNING -> COMPLETED/FAILED）。

## 文档

每个模块和类都提供了详尽的 docstring。请参考代码注释获取更多细节。
