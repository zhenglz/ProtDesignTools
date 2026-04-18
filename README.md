
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
  - `base_tool.py`: 工具基类，定义标准 I/O 范式，支持 CLI 解析和动态传参
- `tools/`: 具体工具实现（如 ProteinMPNN, AF3, Chai-1, ESM2, DLKcat, Pythia, TMalign 等）
- `protocols/`: 流程化协议示例
- `data/`: 配置文件与示例数据

## 快速开始

所有的工具都可以通过以下三种方式进行调用：

### 方式 1: 直接使用工具的 Python 脚本 (CLI 参数)
每个工具内部都封装了 `argparse`，允许你直接从命令行传递参数，而无需编写 JSON 文件。

```bash
python tools/proteinmpnn/tool.py --mode scoring --pdb_path ./example.pdb --mutations "A12G" --exec_mode local
```

### 方式 2: 通过统一的 main.py 和 JSON 配置文件
适用于固定流程和复杂参数组合。

```bash
python main.py proteinmpnn --config data/default_configs/proteinmpnn.json
```

### 方式 3: Python API 调用 (Import)
支持直接实例化类，并在实例化或调用时动态传入参数，非常适合集成到更大的 Python 流程中。

```python
from tools.proteinmpnn.tool import ProteinMPNN

# 初始化时可以指定环境、工作目录等配置
tool = ProteinMPNN(exec_mode="slurm")

# 调用时动态传入运行参数
results = tool(
    mode="design",
    pdb_path="example.pdb",
    num_seqs=5
)
print(results)
```

## 扩展新工具

1. 在 `tools/` 下创建新文件夹（如 `mytool/tool.py`）。
2. 继承 `BaseTool` 类。
3. 重写 `get_cli_parser` 方法，添加该工具特有的命令行参数。
4. 实现 `run` 方法，处理 `input_params` 字典并返回结果字典。
5. 在文件末尾添加 `if __name__ == "__main__": MyTool.cli()`。
6. 在 `main.py` 的字典中注册新工具。

## 任务管理与全局配置

项目包含一个全局配置文件 `data/global_config.json`，控制所有工具的默认行为：
```json
{
    "global_work_dir": "./work_dir",
    "global_exec_mode": "local",
    "default_slurm_params": {
        "partition": "AMD",
        "nodes": 1,
        "ntasks": 1,
        "cpus_per_task": 4
    }
}
```
- **工作目录机制**：每个工具的最终输出目录默认由 `global_work_dir` + `tool_work_dir` (如 `proteinmpnn`) 拼接而成。
- **参数继承**：工具专有的 Slurm 配置会覆盖全局的 `default_slurm_params`，而命令行的临时参数又会覆盖工具的配置。
- **本地模式**：使用 `subprocess` 管理进程。
- **Slurm 模式**：自动生成 `.sh` 脚本，提交并追踪状态（PENDING -> RUNNING -> COMPLETED/FAILED）。

## 文档

每个模块和类都提供了详尽的 docstring。请参考代码注释获取更多细节。
