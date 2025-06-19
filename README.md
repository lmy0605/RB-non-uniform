# `calculateSystemParameters.m` 函数功能说明

## 1. 功能

`calculateSystemParameters.m` 是一个 MATLAB 函数，其核心功能是根据一组基本的物理和计算参数（如瑞利数、网格尺寸等），计算并生成一套用于格子玻尔兹曼方法（Lattice Boltzmann Method, LBM）数值模拟的详细参数。

最终，函数返回一个包含所有计算结果的结构体 `params`，并可以选择性地将这些参数输出到一个日志文件中，以便于记录和追溯。

## 2. 函数接口

```matlab
params = calculateSystemParameters(nx, ny, Rayleigh, Prandtl, constA, logFileName)
```

### 输入参数 (Inputs)

| 参数名           | 类型        | 描述                                                               |
|:------------- |:--------- |:---------------------------------------------------------------- |
| `nx`          | `integer` | x 方向的网格点数。                                                       |
| `ny`          | `integer` | y 方向的网格点数。                                                       |
| `Rayleigh`    | `double`  | 瑞利数 (Ra)，描述浮力驱动与粘性耗散和热耗散之比的关键无量纲数。                               |
| `Prandtl`     | `double`  | 普朗特数 (Pr)，描述动量扩散率与热扩散率之比的无量纲数。                                   |
| `constA`      | `double`  | 一个用于生成非均匀网格的常数，控制网格在边界附近的加密程度。                                   |
| `logFileName` | `string`  | **(可选)** 日志文件的名称。如果提供此参数，所有计算出的参数将被写入该文件；否则，它们将被打印到 MATLAB 命令窗口。 |

### 输出参数 (Output)

| 参数名      | 类型       | 描述                                               |
|:-------- |:-------- |:------------------------------------------------ |
| `params` | `struct` | 一个包含所有计算出的参数的结构体。这个结构体可以方便地传递给其他函数，用于后续的计算和数据处理。 |

好的，这是一个关于您提供的 MATLAB 代码的 Markdown 格式说明文档。

---

# fig-1-2-instT-instU

## `instTv.m`说明

### 1. 功能

读取一个特定时间点的结果二进制文件，对数据进行初步处理（无量纲化），然后将处理后的二维流场数据（速度、温度等）导出为 Tecplot 可视化软件支持的 `.plt` 格式文件。

脚本主要针对二维非均匀网格下的瑞利-伯纳德对流 (Rayleigh-Bénard convection) 问题进行后处理。

### 2. 脚本执行流程

#### 步骤 1: 参数配置 (`basic settings`)

- **文件范围定义**: 设置要处理的文件的起始 (`fileNumStart`)、结束 (`fileNumEnd`) 编号和间隔。
- **路径与命名**: 指定存放二进制数据文件的目录 (`inputDir`) 和文件名的基本部分 (`namebase`)。
- **物理与网格参数**: 定义网格尺寸 (`nx`, `ny`)、拉伸网格常数 (`constA`) 以及流动的物理参数，如瑞利数 (`Rayleigh`) 和普朗特数 (`Prandtl`)。
- **系统参数计算**: 调用一个外部函数 `calculateSystemParameters` 来计算网格坐标 (`xGrid`, `yGrid`) 和单位换算因子 (`velocityUnit`) 等。

#### 步骤 2: 数据读取与处理 (`calculation`)

1. **创建网格坐标**: 使用 `ndgrid` 生成用于 Tecplot 输出的二维坐标矩阵 `Cx` 和 `Cy`。
2. **选取瞬时文件**: 脚本选择指定文件范围中间时刻 (`fileNumStart+ceil(fileSum/2)`) 的一个数据文件作为代表性瞬时场进行处理。
3. **读取二进制数据**: 调用自定义的 `readBinaryFile` 函数，从文件中读取四个物理量的一维数组：
   - `U`: x方向速度
   - `V`: y方向速度
   - `T`: 温度
   - `rho`: 密度
4. **数据重塑与归一化**:
   - 使用 `reshape` 将一维数组转换为 `nx` x `ny` 的二维矩阵。
   - 将速度分量 `U` 和 `V` 除以 `params.velocityUnit`，将其从计算单位转换为无量纲单位。
5. **计算导出量**: 计算速度模 `u_module`，即 `sqrt(U^2 + V^2)`。

#### 步骤 3: 数据导出 (`tec_file`)

- **使用 `liton_ordered_tec` 工具箱**: 这是一个专门用于生成 Tecplot 文件的自定义工具箱。
- **设置导出文件**:
  - 定义输出文件名 `inst_uvT.plt`。
  - 指定要写入文件的变量名：`'X', 'Y', 'U', 'V', 'T', 'Umodule'`。
  - 将处理好的数据（`Cx`, `Cy`, `U`, `V`, `T`, `u_module`）打包。
- **写入文件**: 调用 `.write_plt()` 方法，生成最终的 Tecplot 文件。

### 3. 自定义函数说明

脚本末尾定义了三个辅助函数：

#### a. `readBinaryFile(file, nx, ny)`

- **功能**: 读取特定格式的二进制文件。这种文件格式通常由 Fortran 或 C++ 编写的仿真程序生成，其中数据块之间由整数“标签”隔开。
- **状态**: **此函数在主脚本中被调用**。

#### b. `nonUniformAverage(U, xGrid, yGrid)`

- **功能**: 为在非均匀网格上的场量 `U` 计算基于网格间距的加权平均值。
- **状态**: **注意：此函数在当前脚本中仅被定义，并未被调用**。它可能是为其他分析（如计算全场平均动能）准备的。

#### c. `deri_Twall(T, ...)`

- **功能**: 计算底部和顶部壁面处的温度梯度。这通常用于计算努塞尔数 (Nusselt number)，是热对流问题中的一个关键无量纲参数。该函数使用了适用于拉伸网格的二阶精度的有限差分格式。
- **状态**: **注意：此函数在当前脚本中也仅被定义，并未被调用**。

### 4. 依赖项

要成功运行此脚本，需要确保以下条件满足：

1. **`calculateSystemParameters.m`**: 存在一个名为 `calculateSystemParameters.m` 的函数文件，并且在 MATLAB 的搜索路径中。
2. **`liton_ordered_tec` 工具箱**: 存在一个名为 `liton_ordered_tec` 的类或工具箱，并且在 MATLAB 的搜索路径中。
3. **数据文件**: 在 `inputDir` 指定的路径下，存在格式正确的 `.bin` 二进制数据文件。