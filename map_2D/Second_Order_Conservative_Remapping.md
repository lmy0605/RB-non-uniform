<br />

#### 1. 场重构 (Field Reconstruction)

对于任意一个源单元 `(ii, jj)`，我们不再认为其内部的值处处等于中心值 $Q_{ii,jj}$。相反，我们基于中心值和梯度 $\nabla Q_{ii,jj}$ 来重构一个线性函数 $Q(x,y)$，它近似地描述了该单元内部的物理量分布：

$$
Q(x,y) = Q_{ii,jj} + (x - x_{ii}) \cdot \frac{\partial Q}{\partial x}\bigg|_{ii,jj} + (y - y_{jj}) \cdot \frac{\partial Q}{\partial y}\bigg|_{ii,jj}
$$

其中：

* $Q_{ii,jj}$ 是源单元 `(ii, jj)` 中心的物理量值（例如 `u_src(ii,jj)`）。

* $(x_{ii}, y_{jj})$ 是该源单元的中心坐标。

* $\frac{\partial Q}{\partial x}$ 和 $\frac{\partial Q}{\partial y}$ 是在源单元中心计算出的梯度（即代码中的 `grad_x` 和 `grad_y`）。

#### 2. 守恒积分 (Conservative Integration)

守恒的精髓在于：**一个目标单元** **`(i, j)`** **中物理量的总和，等于所有与之重叠的源单元贡献给它的物理量之和**。

一个源单元 `(ii, jj)` 对一个目标单元 `(i, j)` 的贡献量，就是我们重构的线性函数 $Q(x,y)$ 在它们的**重叠区域** ($A_{\text{overlap}}$) 上的积分：

$$
\text{Contribution} = \iint_{A_{\text{overlap}}} Q(x,y) \,dA
$$

对于一个线性函数在一个矩形区域上的积分，有简化：积分值等于**该区域几何中心点的值**乘以**该区域的面积**。

$$
\iint_{A_{\text{overlap}}} Q(x,y) \,dA = Q(x_c, y_c) \cdot A_{\text{overlap}}
$$

其中 $(x_c, y_c)$ 是重叠区域 $A_{\text{overlap}}$ 的几何中心坐标。代码中 `u_recon * overlap_area` 这行计算的数学基础。

#### 4. 求和与归一化 (Summation & Normalization)

目标单元 `(i, j)` 的最终值 $Q_{\text{tgt}, i,j}$ 是所有源单元贡献的总和，再除以目标单元的总面积 $A_{\text{tgt}, i,j}$ 得到平均值。

$$
Q_{\text{tgt}, i,j} = \frac{\sum_{ii,jj} \left( \iint_{A_{\text{overlap}}(i,j; ii,jj)} Q_{ii,jj}(x,y) \,dA \right)}{A_{\text{tgt}, i,j}}
$$

将第3步的简化代入，得到最终的计算公式：

$$
Q_{\text{tgt}, i,j} = \frac{\sum_{ii,jj} \left( Q_{ii,jj}(x_c, y_c) \cdot A_{\text{overlap}}(i,j; ii,jj) \right)}{\sum_{ii,jj} A_{\text{overlap}}(i,j; ii,jj)}
$$

***

