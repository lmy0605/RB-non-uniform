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

---

# fig5-etaK

### `etaK_timeAvg.m` 说明

分别计算了：

1. **总动能耗散率** (`EUavg`)：基于瞬时速度场 `u` 计算耗散，然后进行时间平均。
2. **湍动能耗散率** (`dissipationAvg`)：基于脉动速度场 `u'` 计算耗散，然后进行时间平均。

### 使用瞬时速度场计算 (Kinetic Energy Dissipation)

这种方法计算的是**总动能的瞬时耗散率 `ε_u`**，然后对其进行时间平均得到 `<ε_u>`。

#### 1. 公式说明

![](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig5formula1.png?raw=true)

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

这种方法计算的是**湍动能 (Turbulent Kinetic Energy, TKE) 的耗散率 `ε`**。

#### 1. 公式说明

![](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig5formula2.png?raw=true)

#### 2. 代码实现

在第二个循环中：

1. `UPrime=U-UAvg; VPrime=V-VAvg;` 计算了脉动速度。
2. `[UX_prime,UY_prime,VX_prime,VY_prime]=GRAD1(UPrime,VPrime,params.dx,params.dy);` 计算了脉动速度的梯度。
3. `dissipation=dissipation+(UX_prime.^2+VY_prime.^2+0.5*(VX_prime+UY_prime).^2);`
   - 这一步计算了括号内的项 `(...)`，并进行了累加。
4. `dissipation=dissipation/fileSum;`
   - 这一步完成了时间平均 `<...>`。
5. `dissipationAvg=viscosity*dissipation*2;`
   - 最后乘以 `2ν`。

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

### `grid_resolution_analysis.m` 说明

### 1. 功能

验证网格尺度小于Kolmogorov尺度。

### 2. 分辨率指标计算

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

---

# fig6-Rascaling

### `ReNu_Ra.m`**说明**

计算全局Nu数和Re数随Ra数变化的**标度率**。

对`Re`、`Nu_wall`（壁面Nu）、`Nu_vol`（体积Nu）、`Nu_eu`（动能耗散Nu）、`Nu_et`（热耗散Nu）与`Ra`数的关系，进行**幂律拟合（Power-Law Fitting）**，确定标度律关系式：

$$
Y = C \cdot Ra^{\alpha}
$$

拟合结果如下：

![](https://github.com/lmy0605/RB-non-uniform/blob/master/screenshots/fig6result.png?raw=true)

---

# fig8-PDF

### `PDF_Nu1e9_fitted.m`说明

对Nu数时间序列数据进行统计分析。

主要功能包括：

- **数据加载**: 从指定的 `.plt` 文件中读取多种Nusselt数（`NuWallAvg`, `NuVolAvg`, `NuEUAvg`, `NuETAvg`）的时间序列。
- **基本统计**: 计算并记录每个Nusselt数时间序列的均值和标准差。
- **概率密度函数 (PDF) 计算**: 对标准化后的数据计算经验PDF（Empirical PDF）。
- **分布拟合**: 将每个Nusselt数脉动数据分别与高斯（Gaussian）分布和广义极值（GEV）分布进行拟合。
- **模型比较**: 通过可视化的PDF曲线对比和量化的信息准则（AIC/BIC），判断哪种分布能更好地描述数据。
- **结果输出**:
  - 将所有统计和分析结果保存到一个详细的日志文件 (`.txt`)。
  - 将PDF数据和拟合曲线数据导出为Tecplot格式的 `.plt` 文件。

---

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

使用 liton_ordered_tec 类将存储了所有 LSC 中心坐标的数组写入 lsc_trajectory.plt 文件。
