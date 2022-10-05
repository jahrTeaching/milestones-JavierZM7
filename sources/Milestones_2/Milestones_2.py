from Modulos.Cauchy_module import Cauchy
from Modulos.Functions import Kepler
from Modulos.Scheme_module import RK4,Euler,Inverse_Euler,Crank_Nicolson

from numpy import array, linspace
import matplotlib.pyplot as plt

U0 = array([1.,0.,0.,1.])
linspace_t = linspace(0,100,10001)                                          # linspace(i,f,d) el salto se calcula como (f-i)/(d-1) #N = linspace_t.shape[0]-1   # los pasos son las particiones de t -1

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

## si se aumentan las particiones en t, hasta los esquemas mas imprecisos se vuelven fiables
# mismas condiciones iniciales
linspace_t = linspace(0,100,100001) # linspace(i,f,d) el salto se calcula como (f-i)/(d-1)

#N = linspace_t.shape[0]-1        #los pasos son las particiones de t -1

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
linspace_t = linspace(0,100,201) # linspace(i,f,d) el salto se calcula como (f-i)/(d-1)

#N = linspace_t.shape[0]-1        #los pasos son las particiones de t -1

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
