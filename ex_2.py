""" 这里尝试使用python搭建一个计算一维热传导方程的求解器
$\frac{\partial u}{\partial t}-a\frac{\partial^2 u}{\partial t^2}=f(x,t),\forall x\in (0,1),t>0 $
$u(x,0)=\varphi(x),\forall x in [0,1]$
$u(0,t)=\alpha(t),u(1,t)=\beta(t),\forall t>0$
其中$f(x,t)=,\varphi(x)=,\alpha(t)=,\beta(t)=$
为此我们构造其差分格式计算T时间内的温度变化情况，记$\Omega=|\{(x,t)| 0<x<1,0<t<T\}$
Step1：使用空间步长$\Delta x==\frac{1}{N}$，时间步长为$\Delta t=\frac{T}{M}$的网格线划分$\Omega$得到网格$\mathbb{J}=\{(x_i,t_j)|(x_i,t_j)\in\bar{\Omega}\}$
    $\mathbb{J}$上$u,f,\varphi,\alpha,\beta$的取值分别记为$u_{ij},f_{ij},\varphi_{ij},\alpha_{ij},\beta_{ij}$
Step2:采用中心差分构造差分格式：
    $\frac{\partial^2 u}{\partial x^2}|_{(x_i,t_j)}=\frac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^2}$
    $\frac{\partial u}{\partial t}|_{x_i,t_j}=\frac{u_{i,j+1}-u_{i,j}}{\Delta t}$
由此得到热传导方程的差分格式：
$\frac{T}{M}(u_{i,j+1}-u_{i,j})-\frac{a}{h^2}(u_{i-1,j}-2u_{i,j}+u_{i+1,j}=f_{i,j},\forall (x_i,t_j)\in \Omega\backlash\partiall\Omega$
边界条件：$u_{i,0}=\varphi_{i},u_{0,j}=\alpha_{j},u_{1,j}=\beta_{j}$
Step3：求解微分方程： 
    $u_{i,j+1}=ru_{i-1,j}+(1-2r)u_{i,j}-ru_{i+1,j}+f_{i,j}, \forall (x_i,t_j)\in \Omega$
    $u_{1,j+1}=(1-2r)u_{1,j}-ru_{2,j}+f_{1,j}+r\alpha_{j}$
    $u_{N-1,j+1}=(1-2r)u_{N-1,j}-ru_{N-2,j}+f_{N-1,j}+r\beta_{j}$  
    $u_{i,0}=\varphi_{i}$"""




import numpy as np
import numpy.matlib
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def var(x):
    #return cos(pi*x)**2*cos(pi*y)**2
    return e**x
def f(x,t):
    #return sin(pi*x)*sin(pi*y)
    return 0
def alp(t):
    return e**t
def beta(t):
    return e**(1+t)
def get_f_vector(j,N,h,s):
    vector=np.matlib.zeros((N-1,1), dtype=float)
    for i in range(0,N-1):
        vector[i]=s*f((i+1)*h,j*s)
    return vector


l=1 #x轴长度
N=20 #x轴等分
T=1 #时间秒
M=2000 #时间等分
a=1 #方程中的系数,热传导系数
s=T/M
h=l/N
r=a*s/(h**2) # r必须小于0.5 此差分方程才稳定
print(r)
x=np.linspace(0,l,N+1, endpoint=True,dtype=float)#
t=np.linspace(0,T,M+1, endpoint=True,dtype=float)
u=np.matlib.zeros((N+1,M+1), dtype=float)
A=np.matlib.identity(N-1, dtype=float)
A=(1-2*r)*A
for i in range(0,N-2):
    A[i,i+1]=r
    A[i+1,i]=r
for i in range(0,N+1):
    u[i,0]=var(i*h)
for j in range(0,M+1):
    u[0,j]=alp(j*s)
    u[N,j]=beta(j*s)
for k in range(1,M+1):
    u[1:N,k]=A*u[1:N,k-1]+get_f_vector(k,N,h,s)
    u[1,k]=u[1,k]+r*u[0,k-1]
    u[N-1,k]=u[N-1,k]+r*u[N,k-1]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
X, Y = np.meshgrid(x, t)
ax.plot_surface(X, Y,u.T , rstride=1, cstride=1, cmap='hot')
plt.show() 

