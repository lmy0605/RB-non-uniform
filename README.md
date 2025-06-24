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

---

# fig-1-2-instT-instU

## `instTv.m`说明

### 1. 功能

读取一个特定时间点的结果二进制文件，对数据进行无量纲化，导出为 Tecplot 可视化软件支持的 `.plt` 格式文件。

### 2. 流程

#### 步骤 1: 参数配置 (`basic settings`)

#### 步骤 2: 数据读取与处理 (`calculation`)

调用自定义的 `readBinaryFile` 函数。将速度分量 `U` 和 `V` 除以 `params.velocityUnit`，将其从计算单位转换为无量纲单位。

#### 步骤 3: 数据导出 (`tec_file`)

- **使用 `liton_ordered_tec` 工具箱**: 这是一个专门用于生成 Tecplot 文件的自定义工具箱。

### 3. 自定义函数说明

#### a. `readBinaryFile(file, nx, ny)`

- **功能**: 读取特定格式的二进制文件。

#### b. `nonUniformAverage(U, xGrid, yGrid)`

- **功能**: 为在非均匀网格上的场量 `U` 计算基于网格间距的加权平均值。
- **状态**: **此函数在当前脚本中仅被定义，并未被调用**。

#### c. `deri_Twall(T, ...)`

- **功能**: 计算底部和顶部壁面处的温度梯度。
- **状态**: **此函数在当前脚本中也仅被定义，并未被调用**。

### 4. 依赖项

1. **`calculateSystemParameters.m`**: 存在一个名为 `calculateSystemParameters.m` 的函数文件，并且在 MATLAB 的搜索路径中。
2. **`liton_ordered_tec` 工具箱**: 存在一个名为 `liton_ordered_tec` 的类或工具箱，并且在 MATLAB 的搜索路径中。
3. **数据文件**: 在 `inputDir` 指定的路径下，存在格式正确的 `.bin` 二进制数据文件。

---

# fig3-timeSeries0

## `timeseries0_2bin.m` 说明

### 1. 功能

计算并导出**雷诺数 (Re)** 和**努塞尔数 (Nu)** 的时间演化。

### 2. 流程

#### 步骤 1: 参数配置

#### 步骤 2: 目录验证

- 开始计算之前，检查指定的数据输入目录是否存在。

#### 步骤 3: 循环与计算

1. **分批次循环处理**
2. **循环内**: 
   - 根据当前数据集的时间步长和文件编号，计算并记录当前的**物理时间**。
   - 调用 `readBinaryFile` 函数读取瞬时速度场 (`U`, `V`) 和温度场 (`T`)。
   - 调用 `deri_Twall` 计算壁面温度梯度。
   - 调用 `nonUniformAverage` 对速度平方和壁面温度梯度进行空间加权平均。

#### 步骤 4: 后处理与导出

1. **累积统计计算**: 调用 `cumulativeAverage` 和 `cumulativePopulationVariance` 函数计算 Re 和 Nu 的**累积平均值**和**累积方差**。
2. **使用 `liton_ordered_tec` 工具箱导出**

### 3. 自定义函数

#### a. `readBinaryFile(file, nx, ny)`

- **功能**: 读取特定格式二进制文件。

#### b. `nonUniformAverage(U, xGrid, yGrid)`

- **功能**: 在非均匀网格上，计算场量 `U` 乘以单元面积（或长度）后的加权值。

#### c. `deri_Twall(T, ...)`

- **功能**: 使用适用于拉伸网格的二阶精度有限差分格式，计算壁面处的温度梯度。

#### d. `cumulativeAverage(U)`

- **功能**: 计算一个向量的累积平均值。

#### e. `cumulativePopulationVariance(U)`

- **功能**: 计算向量的累积总体方差。

### 4. 依赖项

1. **`calculateSystemParameters.m`**: 存在一个名为 `calculateSystemParameters.m` 的函数文件，并且该文件位于 MATLAB 的搜索路径中。
2. **`liton_ordered_tec` 工具箱**: 存在一个名为 `liton_ordered_tec` 的类或工具箱，并且它位于 MATLAB 的搜索路径中，用于生成 `.plt` 文件。

---

# fig4-Nu-timeSeries-stationary

### `timeseries_stationary.m`  说明

### 1. 功能

计算**努塞尔数 (Nu)**

1. ${Nu}_{vol} = \sqrt{{Ra}{Pr}} {\langle u^* T^* \rangle}_V + 1$
2. $Nu_{kinetic}=\sqrt{{Ra}{Pr}} \left\langle \varepsilon_{u}^* \right\rangle_{V}+1$
3. $Nu_{{thermal}} = \sqrt{RaPr} \langle \varepsilon ^{*} _{T}  \rangle_{V}$

脚本生成Nu值的瞬时时间序列（自统计稳态后）及其累积统计特性。

### 2. 代码说明

![](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig4formula.png?raw=true )

已知：

$\left\langle \varepsilon_{\mathfrak{u}} \right\rangle_{V,t} = \frac{\nu^3}{H^4} (Nu - 1) Ra Pr^{-2}$

$\varepsilon_{\mathfrak{u}}(\mathbf{x}, t) = \frac{1}{2} \nu \sum_{ij} \left[ \frac{\partial u_j(\mathbf{x}, t)}{\partial x_i} + \frac{\partial u_i(\mathbf{x}, t)}{\partial x_j} \right]^2$

定义无量纲方式：

$\nu ^*=\sqrt{Pr/Ra}$

$\varepsilon_{\mathfrak{u}}^* (\mathbf{x}, t) = \frac{1}{2} \nu^* \sum_{ij} \left[ \frac{\partial u^*_j(\mathbf{x}, t)}{\partial x^*_i} + \frac{\partial u^*_i(\mathbf{x}, t)}{\partial x^*_j} \right]^2$

可推得$Nu_{kinetic}=\sqrt{{Ra}{Pr}} \left\langle \varepsilon_{u}^* \right\rangle_{V}+1$.

```matlab
[UX,UY,VX,VY]=GRAD1(U,V,params.dx,params.dy);
EUsum=(2.*UX).^2+(2.*VY).^2+2.*(VX+UY).^2;

[~,~,NuEU]=nonUniformAverage(EUsum,params.xGrid,params.yGrid);

NuEUAvg(1,t)=sum(NuEU(:))/params.length0.^2*0.5*params.viscosity0*sqrt(3)/0.1/params.length0*sqrt(Prandtl*Rayleigh)+1;
```

`params.viscosity0`无量纲方式为：$\nu_0=Ma*l_0*\sqrt{Pr/Ra/3}$

即$\nu^*=\nu_0*\sqrt{3}/Ma/l_0$.

所以代码中$\varepsilon^*_{\mathfrak{u}}(\mathbf{x}, t)$计算公式为：

$\varepsilon^*_{\mathfrak{u}}(\mathbf{x}, t) = \frac{1}{2} \nu_0\sqrt{3}/Ma/l_0 \sum_{ij} \left[ \frac{\partial u^*_j(\mathbf{x}, t)}{\partial x^*_i} + \frac{\partial u^*_i(\mathbf{x}, t)}{\partial x^*_j} \right]^2$.

```matlab
%NuEUAvg(1,t)=sum(NuEU(:))/params.length0.^2*0.5*Prandtl+1;
```

或者消去$\nu^*$得到

$Nu_{kinetic}=0.5Pr\sum_{ij} \left[ \frac{\partial u^*_j(\mathbf{x}, t)}{\partial x^*_i} + \frac{\partial u^*_i(\mathbf{x}, t)}{\partial x^*_j} \right]^2$.$

### 3. 自定义函数

#### a. `readBinaryFile(file, nx, ny)`

读取特定格式的二进制文件

#### b. `nonUniformAverage(U, xGrid, yGrid)`

在非均匀网格上，对场量 `U` 进行体积（面积）加权积分的准备步骤。

#### c. `cumulativeAverage(U)`

计算向量的累积平均值。

#### d. `cumulativePopulationVariance(U)`

计算向量的累G积总体方差。

#### e. `GRAD1(U,V,dx,dy)`

在非均匀网格上，使用二阶精度的中心、向前和向后差分格式计算速度场和温度场的一阶梯度。

### 4. 依赖项

1. **`calculateSystemParameters.m`**: 必须存在此函数文件，并位于 MATLAB 搜索路径中。
2. **`liton_ordered_tec` 工具箱**: 必须存在此工具箱，并位于 MATLAB 搜索路径中，用于生成 `.plt` 文件。
