from numpy import linspace,array
import matplotlib.pyplot as plt
from Modulos.Functions import Kepler
from Modulos.Scheme_module import RK4, Crank_Nicolson, Euler, Inverse_Euler
from Error_Richardson import Convergence_R, Error_R

linspace_t = linspace(0,20,2000) #linspace_t = linspace(0,20,2000)
U0 = array([1,0,0,1])

E = Error_R(Kepler,U0,linspace_t,Euler)
E1 = Error_R(Kepler,U0,linspace_t,Inverse_Euler)
E2 = Error_R(Kepler,U0,linspace_t,Crank_Nicolson)
E3 = Error_R(Kepler,U0,linspace_t,RK4)

Ex = E[0,:]
Ey = E[1,:]
Ex1 = E1[0,:]
Ey1 = E1[1,:]
Ex2 = E2[0,:]
Ey2 = E2[1,:]
Ex3 = E3[0,:]
Ey3 = E3[1,:]

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(linspace_t, Ex)
axs[0, 0].plot(linspace_t, Ey)
axs[0, 0].plot(linspace_t, (Ey**2+Ex**2)**0.5)
axs[0, 0].set_title('Euler')
axs[0, 1].plot(linspace_t, Ex1)
axs[0, 1].plot(linspace_t, Ey1)
axs[0, 1].plot(linspace_t, (Ey1**2+Ex1**2)**0.5)
axs[0, 1].set_title('Inverse Euler')
axs[1, 0].plot(linspace_t, Ex2)
axs[1, 0].plot(linspace_t, Ey2)
axs[1, 0].plot(linspace_t, (Ey2**2+Ex2**2)**0.5)
axs[1, 0].set_title('Crank Nicolson')
axs[1, 1].plot(linspace_t, Ex3)
axs[1, 1].plot(linspace_t, Ey3)
axs[1, 1].plot(linspace_t, (Ey3**2+Ex3**2)**0.5)
axs[1, 1].set_title('Runge Kutta 4')

plt.show()

g1 = Convergence_R(Kepler, U0, linspace_t, Euler) #Kepler,u0,lins,Euler
""" print(g1[2])
plt.plot(g1[1],g1[0])
plt.show() """

g2 = Convergence_R(Kepler, U0, linspace_t, Inverse_Euler)
""" print(g2[2])
plt.plot(g2[1],g2[0])
plt.show() """

g3 = Convergence_R(Kepler, U0, linspace_t, Crank_Nicolson)
""" print(g3[2])
plt.plot(g3[1],g3[0])
plt.show() """

g4 = Convergence_R(Kepler, U0, linspace_t, RK4)
""" print(g4[2])
plt.plot(g4[1],g4[0])
plt.show() """

print('q de Euler, Euler_inv ,CN y RK4 ',[g1[2],g2[2],g3[2],g4[2]])

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(g1[1],g1[0])
axs[0, 0].set_title('Euler')
axs[0, 1].plot(g2[1],g2[0])
axs[0, 1].set_title('Inverse Euler')
axs[1, 0].plot(g3[1],g3[0])
axs[1, 0].set_title('Crank Nicolson')
axs[1, 1].plot(g4[1],g4[0])
axs[1, 1].set_title('Runge Kutta 4')
plt.show()