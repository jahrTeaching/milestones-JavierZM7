from numpy import linspace, array, zeros
import matplotlib.pyplot as plt
from Modulos.Cauchy_module import Cauchy
from Modulos.Scheme_module import RK4
from Modulos.Cauchy_module import Cauchy
from Modulos.Functions import KeplerNbody


#introducción de condiciones iniciales #m = 1, G = 1
Nb = 5
#Body1
r01 = array([0.4,0.1,0])
v01 = array([0,0.4,1])
#Body2
r02 = array([0.2,1,0])
v02 = array([0,0.9,1])
#Body3
r03 = array([0.5,0.5,1])
v03 = array([1,-1,0.6])
#Body4
r04 = array([0,0,0])
v04 = array([-0.3,0,-0.3])
#Body5
r05 = array([0.2,-0.1,0.9])
v05 = array([-0.5,1,-0.1])

r0 = zeros([Nb,3]) 
v0 = zeros([Nb,3])

r0[0,:] = r01
r0[1,:] = r02
r0[2,:] = r03
r0[3,:] = r04
r0[4,:] = r05

v0[0,:] = v01
v0[1,:] = v02
v0[2,:] = v03
v0[3,:] = v04 
v0[4,:] = v05

# Ensamblase de matriz inicial
U0 = zeros((2*Nb,3))

for i in range(0,Nb):
    U0[i,:] = r0[i,:]
    U0[Nb+i,:] = v0[i,:]

# Precisión y propagación con RK4

dt = 0.00001
Nt = 200000 
ti = 0
tf = dt*(Nt-1)+ti
lint = linspace(ti,tf,Nt+1)

r = Cauchy(KeplerNbody,U0,lint,RK4)

#Representación de posiciones

ax = plt.figure().add_subplot(projection='3d')
for i in range(0,Nb):
    ax.plot(r[i,0,:],r[i,1,:],r[i,2,:])

ax.set_xlim([-0.5, 1])
ax.set_ylim([-0.2, 1.2])

plt.show() 

