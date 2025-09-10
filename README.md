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

***

# plot\_matlab

### `plot_matlab.m`说明

#### 1. 功能概述

`plot_matlab` 是一个MATLAB绘图函数，可以简化从各种数据源创建图表的过程。

#### 2. 基本用法

函数调用有两种主要形式：

1. **从变量绘图**：
   
   ```matlab
   plot_matlab(x_data, y_data, ...);
   ```

2. **从文件绘图**：
   
   ```matlab
   plot_matlab('Filename', 'path/to/file.plt', ...);
   ```

#### 3. 主要可配置参数详解

以下是 `plot_matlab` 支持的所有参数，可通过 `'ParameterName', ParameterValue` 的形式进行设置。

##### 3.1 数据源参数 (`DataSource` and File-related)

| 参数名              | 类型                          | 默认值          | 描述                                                                                                                                                      |
| ---------------- | --------------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`DataSource`** | `string`                    | `'variable'` | 指定数据来源。可选值为 `'variable'`（从MATLAB变量）或 `'file'`（从文件）。如果指定了 `Filename` 参数，此项会自动设为 `'file'`。                                                                |
| **`Filename`**   | `string` or `cell`          | `''`         | 要读取的数据文件名。可以是单个文件名（字符串），也可以是多个文件名的元胞数组（`{'file1.dat', 'file2.dat'}`）。                                                                                   |
| **`XVariable`**  | `numeric`, `string`, `cell` | `1`          | **仅用于文件模式**。指定用作X轴的数据。可以是：- **数字索引**: `1` (第一列)。- **变量名**: `'Time'`。- **数学表达式**: `'sqrt(X.^2 + Y.^2)'`。- **元胞数组**: `{'Time', 2}` (为多个文件分别指定)。           |
| **`YVariable`**  | `numeric`, `string`, `cell` | *all others* | **仅用于文件模式**。指定用作Y轴的数据。可以是：- **数字索引**: `2` (第二列)。- **变量名**: `'Velocity'`。- **数学表达式**: `'Pressure/1000'`。- **元胞数组**: `{2, 'Pressure', 'U*0.5'}` (绘制多条曲线)。 |

##### 3.2 基础绘图数据 (`x` and `y`)

| 参数名     | 类型                         | 默认值  | 描述                                                                                                                 |
| ------- | -------------------------- | ---- | ------------------------------------------------------------------------------------------------------------------ |
| **`x`** | `vector`                   | `[]` | **仅用于变量模式**。提供X轴的数据向量。                                                                                             |
| **`y`** | `vector`, `matrix`, `cell` | `[]` | **仅用于变量模式**。提供Y轴的数据。可以是：- **向量**: `y` (绘制一条曲线)。- **矩阵**: `[y1; y2]'` (每列是一条曲线)。- **元胞数组**: `{y1, y2}` (每个元胞是一条曲线)。 |

##### 3.3 图例与标签 (`Legend` and `Labels`)

| 参数名                    | 类型        | 默认值     | 描述                                                                  |
| ---------------------- | --------- | ------- | ------------------------------------------------------------------- |
| **`LegendLabels`**     | `cell`    | `{}`    | 一个字符串元胞数组，为每条曲线提供图例标签，例如 `{'Curve 1', 'Curve 2'}`。其数量必须与绘制的曲线数完全匹配。 |
| **`ShowLegend`**       | `logical` | `true`  | 是否显示图例。当有多条曲线时，默认为 `true`。                                          |
| **`XLabel`**           | `string`  | 自动推断    | X轴的标签文本。                                                            |
| **`YLabel`**           | `string`  | 自动推断    | Y轴的标签文本。                                                            |
| **`LabelInterpreter`** | `string`  | `'tex'` | 坐标轴和图例标签的解释器。可选值为 `'tex'`, `'latex'`, `'none'`。                     |

##### 3.4 坐标轴控制 (`Axis Control`)

| 参数名                 | 类型       | 默认值      | 描述                                        |
| ------------------- | -------- | -------- | ----------------------------------------- |
| **`LogScale`**      | `string` | `'none'` | 设置对数坐标轴。可选值为 `'x'`, `'y'`, `'xy'`。        |
| **`XLim`**          | `vector` | `[]`     | 设置X轴的范围，例如 `[0, 10]`。                     |
| **`YLim`**          | `vector` | `[]`     | 设置Y轴的范围，例如 `[-1, 1]`。                     |
| **`XTickInterval`** | `scalar` | `[]`     | **仅用于线性轴**。设置X轴主刻度的间隔。                    |
| **`YTickInterval`** | `scalar` | `[]`     | **仅用于线性轴**。设置Y轴主刻度的间隔。                    |
| **`AxisLayer`**     | `string` | `'top'`  | 将坐标轴层（刻度和框线）置于顶层（`'top'`）或底层（`'bottom'`）。 |
| **`TickDirection`** | `string` | `'out'`  | 设置刻度线的方向，向内（`'in'`）或向外（`'out'`）。          |

##### 3.5 样式与外观 (`Style & Appearance`)

这些参数可以接受**单个值**（应用于所有曲线）或**元胞数组**（循环应用于每条曲线）。

| 参数名                   | 类型                      | 默认值               | 描述                                                                         |
| --------------------- | ----------------------- | ----------------- | -------------------------------------------------------------------------- |
| **`LineColor`**       | `string`, `RGB`, `cell` | 预设颜色              | 线的颜色。支持标准颜色名 (`'r'`), RGB三元组 (`[1 0 0]`), 或预设颜色名 `colorN` (例如 `'color1'`)。 |
| **`LineStyle`**       | `string`, `cell`        | `'-'`             | 线的样式。例如 `'-'`, `'--'`, `':'`, `'-.'`。                                      |
| **`LineWidth`**       | `scalar`, `cell`        | `1.5`             | 线的宽度。                                                                      |
| **`Marker`**          | `string`, `cell`        | `'none'`          | 数据点的标记样式。例如 `'o'`, `'s'`, `'d'`, `'^'`。                                    |
| **`MarkerSize`**      | `scalar`, `cell`        | `8`               | 标记的大小。                                                                     |
| **`MarkerFaceColor`** | `string`, `RGB`, `cell` | `'w'`             | 标记的填充颜色。                                                                   |
| **`MarkerEdgeColor`** | `string`, `RGB`, `cell` | *同* *`LineColor`* | 标记的边缘颜色。                                                                   |
| **`MarkerInterval`**  | `scalar`, `cell`        | `5`               | 每隔多少个数据点显示一个标记。设为`1`显示所有标记。                                                |

##### 3.6 输出参数 (`Output`)

| 参数名                  | 类型       | 默认值                 | 描述                             |
| -------------------- | -------- | ------------------- | ------------------------------ |
| **`OutputFilename`** | `string` | `'plot_output.png'` | 保存输出图像的文件名。                    |
| **`Resolution`**     | `scalar` | `600`               | 输出图像的分辨率 (DPI)。                |
| **`CanvasRatio`**    | `vector` | `[10, 7]`           | 图窗画布的宽高比，例如 `[width, height]`。 |

### `test.m`说明

**示例 1: 从变量绘图，并为不同曲线指定不同样式**

```matlab
x_data = linspace(-5, 5, 200);
y_data = exp(-x_data.^2 / 2) .* cos(3*x_data);
y_data2 = 1.5*cos(3*x_data);
plot_matlab('x', x_data, 'y', {y_data,y_data2}, ...
            'CanvasRatio', [10, 7], ...
            'XLim', [-5, 5], ...
            'YLim', [-2, 2], ...
            'XTickInterval', 2, ...
            'YTickInterval', 0.5, ...
            'LineColor', {[0.8500, 0.3250, 0.0980],...
                     [30/255,144/255,255/255]}, ... % 为每条曲线设置不同颜色
            'LineWidth', 1.5, ...
            'LineStyle', {'-', '-'}, ...  % 为每条曲线设置不同线型
            'Marker', {'o', 's'}, ...   % 为每条曲线设置不同标记
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'w', ...
            'MarkerInterval', 10, ...
            'XLabel', '\it{Time} \rm{(s)}', ...
            'YLabel', '\it{Amplitude}', ...
            'LegendLabels',{'Gauss Modulated','Cosine'},...
            'OutputFilename', 'example1_linear_multi_curve.png', ...
            'Resolution', 600);
```

**示例 2: 从变量绘图，双对数坐标图**

```matlab
x_log = logspace(-1, 2, 100); % 从 10^-1 到 10^2
y_log = 50 * x_log.^-1.5;     % y = 50 * x^(-1.5)
plot_matlab('x', x_log, 'y', y_log, ...
            'LogScale', 'xy', ... % <-- 使用双对数坐标
            'LineWidth', 2, ...
            'Marker', 'o', ...
            'MarkerInterval',10,...
            'MarkerFaceColor', 'w', ...
            'XLabel', '\it{Frequency} \rm{(Hz)}', ...
            'YLabel', 'PSD', ...
            'OutputFilename', 'example2_loglog_plot.png');
```

**示例 3: 从变量绘图，半对数坐标图**

```matlab
x_semilog = linspace(0, 5, 200);
y_semilog = 100 * exp(-1.5 * x_semilog); % 指数衰减
plot_matlab(x_semilog, y_semilog, ...
            'LogScale', 'y', ... % <-- 仅Y轴为对数
            'LineWidth', 2, ...
            'LineColor', 'r', ...
            'XLabel', '\it{t}', ...
            'YLabel', '\it{Concentration} \rm{(log scale)}', ...
            'YLim', [1, 100], ... % 在对数尺度上设置范围
            'OutputFilename', 'example3_semilogy_plot.png');
```

**示例 4: 从`.plt`文件绘图**

```matlab
% 假设 'u_y.plt' 是一个简单的ASCII文件，第一列是时间，其他是数据

if exist('u_y.plt', 'file')
    plot_matlab('DataSource', 'file', 'Filename', 'u_y.plt', ...
                'XVariable', 2, 'YVariable', 3, ...
                'LineColor', 'r', 'XLabel', '\it{y}', 'YLabel', '\it{u^+}','LogScale', 'x', ...
                'OutputFilename', 'example4_file_plot.png');
else
    disp('跳过文件绘图示例: u_y.plt 不存在。');
end
```

**示例 5: 从`.plt`文件绘图，多文件多曲线**

```matlab
  plot_matlab('Filename', {'timeseries_Nu1e8.plt', 'timeseries_Nu1e9.plt','timeseries_Nu1e10.plt','timeseries_Nu1e11.plt','timeseries_Nu1e12.plt'}, ...
              'XVariable', 't', ...
              'XLim',[-0.1,8000],...
              'LineWidth',1,...
              'YVariable', 'NuETAvg', ...
              'XLabel','\it{t}',...
              'YLabel','\it{Nu_{thermal}}',...
              'ShowLegend', 0, ...
              'LogScale','y',...
              'OutputFilename', 'multi_file_comparison.png');
```

**示例 6: 简单计算（+-\*/）**

```matlab
plot_matlab('Filename', 'u_y.plt', ...
            'XVariable', 'Var1', ... % 使用第一列作为X轴
            'YVariable', '(Var3 + Var3+Var4) ', ... % 使用计算出的平均值作为Y轴
            'XLabel', '\it{t}', ...
            'YLabel', 'Sum', ...
            'LegendLabels', {'Sum'}, ... 
            'OutputFilename', 'calculation_from_numeric_file.png');
```

***

# fig1-instT-instU

## `instTv.m`说明

读取一个特定时间点的结果二进制文件，对数据进行无量纲化，导出为 Tecplot 可视化软件支持的 `.plt` 格式文件。

### 2. 流程

#### 步骤 1: 参数配置 (`basic settings`)

在这一步中检查inputDir是否存在。

#### 步骤 2: 数据读取与处理 (`calculation`)

调用自定义的 `readBinaryFile` 函数。将速度分量 `U` 和 `V` 除以 `params.velocityUnit`，将其从计算单位转换为无量纲单位。

#### 步骤 3: 数据导出 (`tec_file`)

**使用** **`liton_ordered_tec`** **工具箱**: 这是一个专门用于生成 Tecplot 文件的自定义工具箱。

### 3. 自定义函数说明

#### a. `readBinaryFile(file, nx, ny)`

* **功能**: 读取特定格式的二进制文件。

#### b. `nonUniformAverage(U, xGrid, yGrid)`

* **功能**: 为在非均匀网格上的场量 `U` 计算基于网格间距的加权平均值。

* **状态**: **此函数在当前脚本中仅被定义，并未被调用**。

#### c. `deri_Twall(T, ...)`

* **功能**: 计算底部和顶部壁面处的温度梯度。

* **状态**: **此函数在当前脚本中也仅被定义，并未被调用**。

### 4. 依赖项

1. **`calculateSystemParameters.m`**: 存在一个名为 `calculateSystemParameters.m` 的函数文件，并且在 MATLAB 的搜索路径中。
2. **`liton_ordered_tec`** **工具箱**: 存在一个名为 `liton_ordered_tec` 的类或工具箱，并且在 MATLAB 的搜索路径中。
3. **数据文件**: 在 `inputDir` 指定的路径下，存在格式正确的 `.bin` 二进制数据文件。

***

# fig3-timeSeries0

## `timeseries0_2bin.m` 说明

计算并导出**雷诺数 (Re)** 和**努塞尔数 (Nu)** 的时间演化。

###### 自定义函数

#### a. `readBinaryFile(file, nx, ny)`

* **功能**: 读取特定格式二进制文件。

#### b. `nonUniformAverage(U, xGrid, yGrid)`

* **功能**: 在非均匀网格上，计算场量 `U` 乘以单元面积（或长度）后的加权值。

#### c. `deri_Twall(T, ...)`

* **功能**: 使用适用于拉伸网格的二阶精度有限差分格式，计算壁面处的温度梯度。

#### d. `cumulativeAverage(U)`

* **功能**: 计算一个向量的累积平均值。

#### e. `cumulativePopulationVariance(U)`

* **功能**: 计算向量的累积总体方差。

***

# fig4-Nu-timeSeries-stationary

### `timeseries_stationary.m`  说明

### 1. 功能

计算**努塞尔数 (Nu)**

1. ${Nu}_{vol} = \sqrt{{Ra}{Pr}} {\langle u^* T^* \rangle}_V + 1$
2. $Nu_{kinetic}=\sqrt{{Ra}{Pr}} \left\langle \varepsilon_{u}^* \right\rangle_{V}+1$
3. $Nu_{{thermal}} = \sqrt{RaPr} \langle \varepsilon ^{*} _{T}  \rangle_{V}$

脚本生成Nu值的瞬时时间序列（自统计稳态后）及其累积统计特性。

### 2. 自定义函数

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

***

# fig5-etaK

### `etaK_timeAvg.m` 说明

分别计算了：

1. **总动能耗散率** (`EUavg`)：基于瞬时速度场 `u` 计算耗散，然后进行时间平均。
2. **湍动能耗散率** (`dissipationAvg`)：基于脉动速度场 `u'` 计算耗散，然后进行时间平均。

### 使用瞬时速度场计算 (Kinetic Energy Dissipation)

这种方法计算的是**总动能的瞬时耗散率** **`ε_u`**，然后对其进行时间平均得到 `<ε_u>`。

#### 1. 公式说明

![1.00](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig5formula1.png?raw=true)

#### 2. 代码实现

在第一个循环中计算 `EUavg`：

```matlab
[UX,UY,VX,VY]=GRAD1(U,V,params.dx,params.dy);
EUavg=EUavg+((2.*UX).^2+(2.*VY).^2+2.*(VX+UY).^2)*viscosity*0.5;
```

Kolmogorov尺度的计算：

$$
\eta_u = \left( \frac{\nu^3}{\langle \epsilon_u \rangle} \right)^{1/4}
$$

代码是 `etaUAvg=(viscosity.^3./EUavg).^0.25;`

### 使用脉动速度场计算 (TKE Dissipation)

这种方法计算的是**湍动能 (Turbulent Kinetic Energy, TKE) 的耗散率** **`ε`**。

#### 1. 公式说明

![1.00](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig5formula2.png?raw=true)

#### 2. 代码实现

在第二个循环中：

1. `UPrime=U-UAvg; VPrime=V-VAvg;` 计算了脉动速度。

2. `[UX_prime,UY_prime,VX_prime,VY_prime]=GRAD1(UPrime,VPrime,params.dx,params.dy);` 计算了脉动速度的梯度。

3. `dissipation=dissipation+(UX_prime.^2+VY_prime.^2+0.5*(VX_prime+UY_prime).^2);`
   
   * 这一步计算了括号内的项 `(...)`，并进行了累加。

4. `dissipation=dissipation/fileSum;`
   
   * 这一步完成了时间平均 `<...>`。

5. `dissipationAvg=viscosity*dissipation*2;`
   
   * 最后乘以 `2ν`。

组合起来，`dissipationAvg` 计算：

$$
2 \nu \left\langle (UX\_prime)^2 + (VY\_prime)^2 + 0.5 \cdot (VX\_prime + UY\_prime)^2 \right\rangle
$$

$$
\epsilon = 2\nu \left\langle \left( \frac{\partial U'}{\partial x} \right)^2 + \left( \frac{\partial V'}{\partial y} \right)^2 + \frac{1}{2} \left( \frac{\partial V'}{\partial x} + \frac{\partial U'}{\partial y} \right)^2 \right\rangle
$$

Kolmogorov尺度的计算：

$$
\eta_k = \left( \frac{\nu^3}{\epsilon} \right)^{1/4}
$$

代码是 `etaKAvg=(viscosity.^3./dissipationAvg).^0.25;`

#### 3. grid\_resolution\_analysis

验证网格尺度小于Kolmogorov尺度。

脚本计算了两种分辨率指标 (`grid_resolution` 和 `grid_resolution2`)。

首先计算:

$$
\Delta x_i=\frac{x_{i+1}-x_{i-1}}{2}
$$

$$
\Delta y_j=\frac{y_{j+1}-y_{j-1}}{2}
$$

* **标准一：最长边**
  
  $$
  \Delta_1(i,j) = \max(\Delta x_i, \Delta y_j)
  $$
  
  ```matlab
  grid_resolution = max(deltaX, deltaY); 
  ```

* **标准二：等效面积**
  
  $$
  \Delta_2(i,j) = \sqrt{\Delta x_i \cdot \Delta y_j}
  $$
  
  此处的 `deltaxy` 是由函数直接返回的面积矩阵 `node_area_weights`。
  
  ```matlab
  grid_resolution2 = sqrt(deltaxy);
  ```

将上述两种网格分辨率分别除以输入的物理尺度 `etaK` 和 `etaU`。

$$
\text{Resolution Ratio}(i,j) = \frac{\Delta_2(i,j)}{\eta_K(i,j)}
$$

```matlab
% 使用等效面积法计算与Kolmogorov动能耗散尺度的比率
resolution_ratio_etaK2 = grid_resolution2 ./ etaK;

% 使用等效面积法计算与Batchelor热耗散尺度的比率
resolution_ratio_etaU2 = grid_resolution2 ./ etaU;
```

展弦比是衡量单元拉伸程度的指标，理想值为1（正方形）。

$$
\text{Aspect Ratio}(i,j) = \frac{\max(\Delta x_i, \Delta y_j)}{\min(\Delta x_i, \Delta y_j)}
$$

```matlab
aspect_ratio = max(deltaX, deltaY) ./ min(deltaX, deltaY);
```

***

# fig6-Rascaling

### `ReNu_Ra.m`**说明**

读取全局Nu数和Re数（Nu\_Ra.dat），计算随Ra数变化的**标度率**。

对`Re`、`Nu_wall`（壁面Nu）、`Nu_vol`（体积Nu）、`Nu_eu`（动能耗散Nu）、`Nu_et`（热耗散Nu）与`Ra`数的关系，进行**幂律拟合（Power-Law Fitting）**，确定标度律关系式：

$$
Y = C \cdot Ra^{\alpha}
$$

拟合结果如下：

![1.00](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig6result.png?raw=true)

***

# fig8-PDF

### `PDF_Nu1e9_fitted.m`说明

对Nu数时间序列数据进行统计分析。

主要功能包括：

* **数据加载**: 从指定的 `.plt` 文件中读取多种Nusselt数（`NuWallAvg`, `NuVolAvg`, `NuEUAvg`, `NuETAvg`）的时间序列。

* **基本统计**: 计算并记录每个Nusselt数时间序列的均值和标准差。

* **概率密度函数 (PDF) 计算**: 对标准化后的数据计算经验PDF（Empirical PDF）。

* **分布拟合**: 将每个Nusselt数脉动数据分别与高斯（Gaussian）分布和广义极值（GEV）分布进行拟合。

* **模型比较**: 通过可视化的PDF曲线对比和量化的信息准则（AIC/BIC），判断哪种分布能更好地描述数据。

* **结果输出**:
  
  * 将所有统计和分析结果保存到一个详细的日志文件 (`.txt`)。
  
  * 将PDF数据和拟合曲线数据导出为Tecplot格式的 `.plt` 文件。

***

# fig12-vor

## `vorticity_inden_sort.m`说明

计算$\omega _ z$,阈值离散并确定中心

## `vorticityBuoyancy_inden_sort.m`说明

计算$\nabla \times F$,阈值离散并确定中心

## `enstro_inden_sort.m`说明

计算$\omega * (\nabla \times F)$,阈值离散并确定中心

## `center_tra_new.m`说明

###### 自动判断 LSC 的旋转方向（顺时针或逆时针）：

通过比较中心区域的正负涡量贡献来实现。

```matlab
    pos_mask = vor_z_direction > 0;
    neg_mask = vor_z_direction < 0;

    positive_enstrophy = sum(vor_z_direction(pos_mask).^2 .* area_weights_direction(pos_mask));
    negative_enstrophy = sum(vor_z_direction(neg_mask).^2 .* area_weights_direction(neg_mask));

    is_LSC_clockwise = negative_enstrophy > positive_enstrophy;
```

###### 识别 LSC中心：

1. 根据涡量方向和设定的阈值生成二值化图像。
2. 识别图像中所有的涡旋结构。使用Image processing toolbox中的`bwlabel`函数。
3. 计算**每一个**涡旋的物理面积和质心，筛选出位于计算域中心区域的涡旋。使用Image processing toolbox中的`regionprops`函数。
4. 在中心区域的涡旋中，选择面积最大的一个作为 LSC。

###### 输出结果：

使用 liton\_ordered\_tec 类将存储了所有 LSC 中心坐标的数组写入 lsc\_trajectory.plt 文件。

***

# fig14-tke

### `tke.m`说明

读取一系列瞬时流场文件（`.bin`），计算时间平均统计量，最终计算出TKE产生项和耗散项的空间分布保存为Tecplot可读的 `.plt` 文件。

* **TKE产生项计算**:
  
  * **剪切产生项 (Shear Production)**: `-<u'_i u'_j> * (∂<U_i>/∂x_j)`
  
  * **浮力产生项 (Buoyancy Production)**: `<v'T'>` (假设重力在y方向)

* **TKE耗散项计算**:
  
  * **伪耗散 (Pseudo-dissipation)**: `ν * <(∂u'_i/∂x_j) * (∂u'_i/∂x_j)>`
  
  * **真耗散 (True dissipation)**: `2 * ν * <s'_{ij}s'_{ij}>`，其中 `s'_{ij}` 是脉动应变率张量。

* **非均匀网格处理**: 使用自定义的有限差分格式（`GRAD1`函数）和积分权重（`nonUniformAverage`函数）来处理非均匀网格。

### `tke_ra.m`说明

读取全场体积平均的产生率和耗散率（tke\_ra.dat）。

* **可视化**:
  
  * 调用 `plot_matlab` 函数，分别生成两张PNG图。
  
  * **图1 (`p_ra.png`)**: TKE产生率 vs. 瑞利数。
  
  * **图2 (`d_ra.png`)**: TKE耗散率 vs. 瑞利数。

***

# fig15-Fourier

### `Fourier_map.m`说明

#### 1. 功能

* **数据读取**: 读取指定文件夹内的二进制速度场文件 (`.bin`)。

* **插值：** 非均匀网格插值为均匀网格。插值函数`conservative_remap_2d.m`采用二阶保守重映射，详细见‘map\_2D/Second\_Order\_Conservative\_Remapping.md’说明。

* **傅里叶分解**: 生成二维正弦/余弦傅里叶基函数，并将每个时间点的速度场投影到这些基函数上，以获得傅里叶系数 (`Ax`, `Ay`)。

* **能量计算**: 基于傅里叶系数，计算每个模态的瞬时能量、总能量、能量百分比，以及能量的时间标准差（`E_rms`）。

* **统计分析**: 计算各模态能量的时间平均值、稳定性参数 (`S11`) 等物理量。

* **数据导出**: 存为ASCII 文本文件 (`.dat`)。

#### 2. 输出文件说明

* `Fourier_coeff_Ax_[modeNum].dat`: 第 `modeNum` 个模态的 `Ax` 系数随时间变化的数据。

* `Fourier_coeff_Ay_[modeNum].dat`: 第 `modeNum` 个模态的 `Ay` 系数随时间变化的数据。

* `energyPercentage_[modeNum].dat`: 第 `modeNum` 个模态的能量百分比随时间变化的数据。

* **`Fourier_result_summary.dat`**:这是一个矩阵，每行代表一个模态，每列代表一个统计量：
  
  * **列 1**: `mean(abs(Ax))`
  
  * **列 2**: `mean(abs(Ay))`
  
  * **列 3**: `<E^(m,n)(t)>` (模态能量的时间平均值)
  
  * **列 4**: `<E^(m,n)(t) / E_total(t)>` (能量百分比的时间平均值)
  
  * **列 5**: `<E^(m,n)(t)> / <E_total(t)> * 100`
  
  * **列 6**: `S11 = <E^(m,n)> / E_rms^(m,n)` (稳定性参数)
  
  * **列 7**: `E_rms^(m,n) / <E_total>`

* `energyPercentage.png` (如果 `exportFigure=1`): 显示所有模态能量百分比随时间变化的图表。
