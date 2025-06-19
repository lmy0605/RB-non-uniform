# `calculateSystemParameters.m` 函数功能说明

## 1. 功能

`calculateSystemParameters.m` 是一个 MATLAB 函数，其核心功能是根据一组基本的物理和计算参数（如瑞利数、网格尺寸等），计算并生成一套用于格子玻尔兹曼方法（Lattice Boltzmann Method, LBM）数值模拟的详细参数。

最终，函数返回一个包含所有计算结果的结构体 `params`，并可以选择性地将这些参数输出到一个日志文件中，以便于记录和追溯。

## 2. 函数接口

```matlab
params = calculateSystemParameters(nx, ny, Rayleigh, Prandtl, constA, logFileName)
```

### 输入参数 (Inputs)

| 参数名 | 类型 | 描述 |
| :--- | :--- | :--- |
| `nx` | `integer` | x 方向的网格点数。 |
| `ny` | `integer` | y 方向的网格点数。 |
| `Rayleigh` | `double` | 瑞利数 (Ra)，描述浮力驱动与粘性耗散和热耗散之比的关键无量纲数。 |
| `Prandtl` | `double` | 普朗特数 (Pr)，描述动量扩散率与热扩散率之比的无量纲数。 |
| `constA` | `double` | 一个用于生成非均匀网格的常数，控制网格在边界附近的加密程度。 |
| `logFileName`| `string` | **(可选)** 日志文件的名称。如果提供此参数，所有计算出的参数将被写入该文件；否则，它们将被打印到 MATLAB 命令窗口。 |

### 输出参数 (Output)

| 参数名 | 类型 | 描述 |
| :--- | :--- | :--- |
| `params` | `struct` | 一个包含所有计算出的参数的结构体。这个结构体可以方便地传递给其他函数，用于后续的计算和数据处理。 |
