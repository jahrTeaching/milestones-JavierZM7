##Sin terminar!
import matplotlib.pyplot as plt
from numpy import linspace,size,zeros,array,linalg,shape,log10

from Modulos.Cauchy_module import Cauchy
from Modulos.Functions import Kepler
from Modulos.Scheme_module import Crank_Nicolson,RK4,Euler,Inverse_Euler   #, Euler, Inverse_Euler,RK4
Yaxis = zeros([3])
Xaxis = zeros([3])
#PUNTO 1
N=2
U0 = array([1,0,0,1])
print('el vector U0 es',U0)
t=linspace(0,100,N)
print('el mallado original, t, es:',t)
t2=linspace(0,100,2*N-1)
print('El mallado 2 es t2:',t2)

Un  = Cauchy(Kepler,U0,t,Euler)
U2n = Cauchy(Kepler,U0,t2,Euler)
print('La matriz Un es',Un)
print('La matriz U2n es',U2n)
#q=1
R=zeros([size(U0),1])
Ev=zeros([N+1])
#print(R)
#R[:,1]=((U2n[:,2*1]-Un[:,1])/(1-(1/2**q)))
#print(R)

R[:,0]=((U2n[:,2*N-2]-Un[:,N-1]))#/(1-(1/2**q)))
Rlin=linalg.norm(R)

print(Rlin) #este es el valor de la norma de E para N=4

Yaxis[0] = log10(Rlin)
Xaxis[0] = log10(N)

#punto2

N=200
U0 = array([1,0,0,1])
print('el vector U0 es',U0)
t=linspace(0,100,N)
print('el mallado original, t, es:',t)
t2=linspace(0,100,2*N-1)
print('El mallado 2 es t2:',t2)

Un  = Cauchy(Kepler,U0,t,Euler)
U2n = Cauchy(Kepler,U0,t2,Euler)
print('La matriz Un es',Un)
print('La matriz U2n es',U2n)
#q=1
R=zeros([size(U0),1])
Ev=zeros([N+1])
#print(R)
#R[:,1]=((U2n[:,2*1]-Un[:,1])/(1-(1/2**q)))
#print(R)

R[:,0]=((U2n[:,2*N-2]-Un[:,N-1]))#/(1-(1/2**q)))
Rlin=linalg.norm(R)

print(Rlin) #este es el valor de la norma de E para N=4

Yaxis[1] = log10(Rlin)
Xaxis[1] = log10(N)

#punto3

N=20000
U0 = array([1,0,0,1])
print('el vector U0 es',U0)
t=linspace(0,100,N)
print('el mallado original, t, es:',t)
t2=linspace(0,100,2*N-1)
print('El mallado 2 es t2:',t2)

Un  = Cauchy(Kepler,U0,t,Euler)
U2n = Cauchy(Kepler,U0,t2,Euler)
print('La matriz Un es',Un)
print('La matriz U2n es',U2n)
q=1
R=zeros([size(U0),1])
Ev=zeros([N+1])
#print(R)
#R[:,1]=((U2n[:,2*1]-Un[:,1])/(1-(1/2**q)))
#print(R)

R[:,0]=((U2n[:,2*N-2]-Un[:,N-1]))#/(1-(1/2**q)))
Rlin=linalg.norm(R)

print(Rlin) #este es el valor de la norma de E para N=4

Yaxis[2] = log10(Rlin)
Xaxis[2] = log10(N)

plt.plot(Xaxis, Yaxis,'.')
plt.show()

## en bucle
R=zeros([size(U0),1])
puntos = 25
Yaxis = zeros([puntos-2])
Xaxis = zeros([puntos-2])
N=1
for i in range (2,puntos):
    N=2*N
    U0 = array([1,0,0,1])
    #print('el vector U0 es',U0)
    t=linspace(0,100,N)
    #print('el mallado original, t, es:',t)
    t2=linspace(0,100,2*N)
    #print('El mallado 2 es t2:',t2)

    Un  = Cauchy(Kepler,U0,t,Euler)
    U2n = Cauchy(Kepler,U0,t2,Euler)
    #print('La matriz Un es',Un)
    #print('La matriz U2n es',U2n)
    #q=1
    
    #Ev=zeros([N+1])
    #print(R)
    #R[:,1]=((U2n[:,2*1]-Un[:,1])/(1-(1/2**q)))
    #print(R)

    R[:,0]=((U2n[:,-1]-Un[:,-1]))#/(1-(1/2**q)))
    Rlin=linalg.norm(R)

    #print(Rlin) #este es el valor de la norma de E para N=4
    Yaxis[i-2] = log10(Rlin)
    Xaxis[i-2] = log10(N)

plt.plot(Xaxis, Yaxis,)
plt.show()
