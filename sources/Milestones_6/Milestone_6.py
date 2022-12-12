from numpy import array, linspace, zeros
import matplotlib.pyplot as plt
from random import random
from Modulos.Cauchy_module import Cauchy
from Modulos.Scheme_module import RKEmb
from Modulos.CR3BPF import Lpoints, StabLP
from Modulos.Functions import MoonEarthCR3BP

##Datos
N =1000
tf = 10
t = linspace(0, tf, N) 
mu = 1.2151e-2 #Tierra-Luna
#mu = 3.0039e-7 #Sol-Tierra

# puntos iniciales para iterar y determinar la posición de los puntos de Lagrange
U0LP = array([[0.5, 0, 0, 0],[1.5, 0, 0, 0],[-0.5, 0, 0, 0],[0.5, 0.5, 0, 0],[0.5, -0.5, 0, 0]])
LagPoints = Lpoints(U0LP, 5, mu)
print(LagPoints)

#Estudio la estabilidad de los puntos de Lagrange

U0S = zeros(4)
U0SLP = zeros((5,4))
for i in range(5):
    U0S = zeros(4)
    U0S[0:2] = LagPoints[i,:]
    eingvalues = StabLP(U0S, mu)
    U0SLP[i,:] = eingvalues.real
print(U0SLP)

#Cálculo de orbita para Cada punto de Lagrange, las condiciones iniciales va a ser un punto cercano al punto de Lagrange, sumando un epsilon
ULG = zeros((4,N,5))
for LG in range(5):
    U0 = zeros(4)
    U0[0:2] = LagPoints[LG-1,:]
    eps = 1e-4*random()
    U0 = U0 + eps
    U = Cauchy(MoonEarthCR3BP, U0, t, RKEmb)
    ULG[:,:,LG-1] = U[:,:]

print(ULG)

#Grafico

fig, axs = plt.subplots(2, 3)
axs[0, 0].plot(ULG[0,:,0], ULG[1,:,0],'-',color = "r")
axs[0, 0].plot(ULG[0,:,1], ULG[1,:,1],'-',color = "r")
axs[0, 0].plot(ULG[0,:,2], ULG[1,:,2],'-',color = "r")
axs[0, 0].plot(ULG[0,:,3], ULG[1,:,3],'-',color = "r")
axs[0, 0].plot(ULG[0,:,4], ULG[1,:,4],'-',color = "r")
axs[0, 0].plot(-mu, 0, 'o', color = "b")
axs[0, 0].plot(1-mu, 0, 'o', color = "#808080")
for i in range(5):
    axs[0, 0].plot(LagPoints[i,0], LagPoints[i,1] , 'o', color = "k")
axs[0, 0].set_title("Lagrange point View")
axs[0, 0].locator_params(axis='x', nbins=5)
axs[0, 0].locator_params(axis='y', nbins=5)

axs[0, 1].plot(ULG[0,:,4], ULG[1,:,4],'-',color = "r")
axs[0, 1].plot(LagPoints[4,0], LagPoints[4,1] , 'o', color = "k")
axs[0, 1].set_title("L5 View")
axs[0, 1].locator_params(axis='x', nbins=5)
axs[0, 1].locator_params(axis='y', nbins=5)

axs[0, 2].plot(ULG[0,:,3], ULG[1,:,3],'-',color = "r")
axs[0, 2].plot(LagPoints[3,0], LagPoints[3,1] , 'o', color = "k")
axs[0, 2].set_title("L4 View")
axs[0, 2].locator_params(axis='x', nbins=5)
axs[0, 2].locator_params(axis='y', nbins=5)

axs[1, 0].plot(ULG[0,:,2], ULG[1,:,2],'-',color = "r")
axs[1, 0].plot(LagPoints[2,0], LagPoints[2,1] , 'o', color = "k")
axs[1, 0].set_title("L3 View")
axs[1, 0].locator_params(axis='x', nbins=5)
axs[1, 0].locator_params(axis='y', nbins=5)

axs[1, 1].plot(ULG[0,:,1], ULG[1,:,1],'-',color = "r")
axs[1, 1].plot(LagPoints[1,0], LagPoints[1,1] , 'o', color = "k")
axs[1, 1].set_title("L2 View")
axs[1, 1].locator_params(axis='x', nbins=5)
axs[1, 1].locator_params(axis='y', nbins=5)

axs[1, 2].plot(ULG[0,:,0], ULG[1,:,0],'-',color = "r")
axs[1, 2].plot(LagPoints[0,0], LagPoints[0,1] , 'o', color = "k")
axs[1, 2].set_title("L1 View")
axs[1, 2].locator_params(axis='x', nbins=5)
axs[1, 2].locator_params(axis='y', nbins=5)

plt.show()
