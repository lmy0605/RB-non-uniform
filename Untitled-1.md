已知：

\left\langle \varepsilon_{\mathfrak{u}} \right\rangle_{V,t} = \frac{\nu^3}{H^4} (Nu - 1) Ra Pr^{-2}

\varepsilon_{\mathfrak{u}}(\mathbf{x}, t) = \frac{1}{2} \nu \sum_{ij} \left[ \frac{\partial u_j(\mathbf{x}, t)}{\partial x_i} + \frac{\partial u_i(\mathbf{x}, t)}{\partial x_j} \right]^2

定义无量纲方式：

\nu^*=\sqrt{Pr/Ra}

\varepsilon^*_{\mathfrak{u}}(\mathbf{x}, t) = \frac{1}{2} \nu^* \sum_{ij} \left[ \frac{\partial u^*_j(\mathbf{x}, t)}{\partial x^*_i} + \frac{\partial u^*_i(\mathbf{x}, t)}{\partial x^*_j} \right]^2

可推得Nu_{kinetic}=\sqrt{{Ra}{Pr}} \left\langle \varepsilon_{u}^* \right\rangle_{V}+1.

```matlab
[UX,UY,VX,VY]=GRAD1(U,V,params.dx,params.dy);
EUsum=(2.*UX).^2+(2.*VY).^2+2.*(VX+UY).^2;

[~,~,NuEU]=nonUniformAverage(EUsum,params.xGrid,params.yGrid);

NuEUAvg(1,t)=sum(NuEU(:))/params.length0.^2*0.5*params.viscosity0*sqrt(3)/0.1/params.length0*sqrt(Prandtl*Rayleigh)+1;
```

`params.viscosity0`无量纲方式为：\nu_0=Ma*l_0*\sqrt{Pr/Ra/3}

即\nu^*=\nu_0*\sqrt{3}/Ma/l_0.

所以代码中\varepsilon^*_{\mathfrak{u}}(\mathbf{x}, t)计算公式为：

\varepsilon^*_{\mathfrak{u}}(\mathbf{x}, t) = \frac{1}{2} \nu_0\sqrt{3}/Ma/l_0 \sum_{ij} \left[ \frac{\partial u^*_j(\mathbf{x}, t)}{\partial x^*_i} + \frac{\partial u^*_i(\mathbf{x}, t)}{\partial x^*_j} \right]^2.

```matlab
%NuEUAvg(1,t)=sum(NuEU(:))/params.length0.^2*0.5*Prandtl+1;
```

或者消去\nu^*得到

$Nu
