# HeatConductionEquation
这里尝试使用python搭建一个计算一维热传导方程的差分方法发求解器

$\frac{\partial u}{\partial t}-a\frac{\partial^2 u}{\partial t^2}=f(x,t),\forall x\in (0,1),t>0 $

$u(x,0)=\varphi(x),\forall x in \[0,1\]$

$u(0,t)=\alpha(t),u(1,t)=\beta(t),\forall t>0$

其中

$f(x,t)=,\varphi(x)=,\alpha(t)=,\beta(t)=$

为此我们构造其差分格式计算T时间内的温度变化情况，记

$\Omega=\{(x,t)| 0<x<1,0<t<T\}$

Step1：使用空间步长

$\Delta x==\frac{1}{N}$

，时间步长为

$\Delta t=\frac{T}{M}$

的网格线划分$\Omega$得到网格

$\mathbb{J}=\{(x_i,t_j) |(x_i,t_j) \in \bar{\Omega} \}$

$\mathbb{J}$

上

$u,f,\varphi,\alpha,\beta$

的取值分别记为

$u_{ij},f_{ij},\varphi_{ij},\alpha_{ij},\beta_{ij}$

Step2:采用中心差分构造差分格式：

$\frac{\partial^2 u}{\partial x^2}=\frac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^2}$

$\frac{\partial u}{\partial t}=\frac{u_{i,j+1}-u_{i,j}}{\Delta t}$

由此得到热传导方程的差分格式：

$\frac{T}{M}(u_{i,j+1}-u_{i,j})-\frac{a}{h^2}(u_{i-1,j}-2u_{i,j}+u_{i+1,j}=f_{i,j},\forall (x_i,t_j)\in \Omega$

边界条件：

$u_{i,0}=\varphi_{i},u_{0,j}=\alpha_{j},u_{1,j}=\beta_{j}$

Step3：求解微分方程： 

$u_{i,j+1}=ru_{i-1,j}+(1-2r)u_{i,j}-ru_{i+1,j}+f_{i,j}, \forall (x_i,t_j)\in \Omega$

$u_{1,j+1}=(1-2r)u_{1,j}-ru_{2,j}+f_{1,j}+r\alpha_{j}$

$u_{N-1,j+1}=(1-2r)u_{N-1,j}-ru_{N-2,j}+f_{N-1,j}+r\beta_{j}$  

$u_{i,0}=\varphi_{i}$
