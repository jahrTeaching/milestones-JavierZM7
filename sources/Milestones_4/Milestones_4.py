from Modulos.Cauchy_module import Cauchy
from Modulos.Functions import Oscilator_lin
from Modulos.Scheme_module import RK4,Euler,Inverse_Euler,Crank_Nicolson,Leap_frog 
from Modulos.Regions import Stab_Regions

from numpy import array, linspace, size
import matplotlib.pyplot as plt

U0 = array([1,0]) 
dt = array([0.01, 0.1, 0.5, 0.9, 1])   
tf = 100    

for i in range(0,size(dt)): 
    N = ((100/dt[i])+1) 
    N = round(N)
    linspace_t = linspace(0,tf,(N))  

    y1 = Cauchy(Oscilator_lin,U0,linspace_t,Euler)
    y2 = Cauchy(Oscilator_lin,U0,linspace_t,Inverse_Euler)
    y3 = Cauchy(Oscilator_lin,U0,linspace_t,Crank_Nicolson)
    y4 = Cauchy(Oscilator_lin,U0,linspace_t,RK4)
    y5 = Cauchy(Oscilator_lin,U0,linspace_t,Leap_frog)

    fig, axs = plt.subplots(2, 3)
    axs[0, 0].plot(linspace_t, y1[0,:])
    axs[0, 0].set_title('Euler')
    axs[0, 1].plot(linspace_t, y2[0,:])
    axs[0, 1].set_title('Inverse Euler')
    axs[1, 0].plot(linspace_t, y3[0,:])
    axs[1, 0].set_title('Crank Nicolson')
    axs[1, 1].plot(linspace_t, y4[0,:])
    axs[1, 1].set_title('Runge Kutta 4')
    axs[1, 2].plot(linspace_t, y5[0,:])
    axs[1, 2].set_title('Leap-Frog')
    plt.show()

p = 1000 #precisión de malla de region
ejex = linspace(-5, 5, p)
ejey = linspace(-3.5, 3.5, p)
#Defino los ejes para representar los complejos
stab_Euler = Stab_Regions(ejex,ejey,Euler)
stab_Inverse_Euler = Stab_Regions(ejex,ejey,Inverse_Euler)
stab_CN = Stab_Regions(ejex,ejey,Crank_Nicolson)
stab_RK4 = Stab_Regions(ejex,ejey,RK4)
stab_LP = Stab_Regions(ejex,ejey,Leap_frog)

fig, axs = plt.subplots(2, 2)
axs[0, 0].contourf(ejex, ejey, stab_Euler, levels = [0, 1],  colors=['#C0C0C0'] )
axs[0, 0].contour(ejex, ejey, stab_Euler, levels = [0.5,1,1.5] )
axs[0, 0].plot([0,0],[dt,-dt], 'o', label = "Raíces del oscilador lineal con dt = " +str(dt))
axs[0, 0].set_title('Estabilidad Euler')
axs[0, 1].contourf(ejex, ejey, stab_Inverse_Euler, levels = [0, 1],  colors=['#C0C0C0'] )
axs[0, 1].contour(ejex, ejey, stab_Inverse_Euler, levels = [0.5,1,1.5] )
axs[0, 1].plot([0,0],[dt,-dt], 'o', label = "Raíces del oscilador lineal con dt = " +str(dt))
axs[0, 1].set_title('Estabilidad Inverse Euler')
axs[1, 0].contourf(ejex, ejey, stab_CN, levels = [0, 1],  colors=['#C0C0C0'] )
axs[1, 0].contour(ejex, ejey, stab_CN, levels = [0.5,1,1.5] )
axs[1, 0].plot([0,0],[dt,-dt], 'o', label = "Raíces del oscilador lineal con dt = " +str(dt))
axs[1, 0].set_title('Estabilidad Crank Nicolson')
axs[1, 1].contourf(ejex, ejey, stab_RK4, levels = [0, 1],  colors=['#C0C0C0'] )
axs[1, 1].contour(ejex, ejey, stab_RK4, levels = [0.5,1,1.5] )
axs[1, 1].plot([0,0],[dt,-dt], 'o', label = "Raíces del oscilador lineal con dt = " +str(dt))
axs[1, 1].set_title('Estabilidad Runge Kutta 4')

plt.show()