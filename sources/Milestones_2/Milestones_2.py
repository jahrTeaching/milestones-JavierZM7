from Modulos.Cauchy_module import Cauchy
from Modulos.Functions import Kepler
from Modulos.Scheme_module import RK4,Euler,Inverse_Euler,Crank_Nicolson #para usar newton, cambiar a Modulos.Scheme_module_Newton

from numpy import array, linspace
import matplotlib.pyplot as plt

U0 = array([1.,0.,0.,1.])
linspace_t = linspace(0,100,10001)                                         

Solución = Cauchy(Kepler,U0,linspace_t,Euler)
#plot 
Xorb = Solución[0,:]
Yorb = Solución[1,:] 
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,Inverse_Euler)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,Crank_Nicolson)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,RK4)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

# mismas condiciones iniciales
linspace_t = linspace(0,100,100001) 

Solución = Cauchy(Kepler,U0,linspace_t,Euler)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:] 
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,Inverse_Euler)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,Crank_Nicolson)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,RK4)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

U0 = array([1.,0.,0.,1.])
linspace_t = linspace(0,100,201) 

Solución = Cauchy(Kepler,U0,linspace_t,Crank_Nicolson)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()

Solución = Cauchy(Kepler,U0,linspace_t,RK4)
#plot de la órbita
Xorb = Solución[0,:]
Yorb = Solución[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()
