import numpy as np
import matplotlib.pyplot as plt

#Condiciones iniciales
x0 = 1.
vx0 = 0.
y0 = 0.
vy0 = 1.
U0 = np.array([[x0],[y0],[vx0],[vy0]])
U0.shape = (4,1)
#saltos de tiempo y tiempo final
t0=0
DeltaT= 0.01
tf=10
#Metodo Euler
F0 = np.array([[vx0],[vy0],[-x0/((x0**2+y0**2.)**(3./2.))],[-y0/((x0**2.+y0**2.)**(3./2.))]])
F0.shape = (4,1)
t=t0
#En este vector se almacena los n resultados en columnas
Orbita2D = np.array([[0.],[0.],[0.],[0.]])
Orbita2D.shape = (4,1)
#Defino vectores
U = np.array([[0.],[0.],[0.],[0.]])
U.shape = (4,1)
Upre = np.array([[0.],[0.],[0.],[0.]])
U.shape = (4,1)
Fpre = np.array([[0.],[0.],[0.],[0.]])
U.shape = (4,1)
Upre=U0
Orbita2D=U0
while (t < tf):
    Fpre = np.array([[Upre[2,0]],[Upre[3,0]],[-(Upre[0,0])/(((Upre[0,0])**2+(Upre[1,0])**2.)**(3./2.))],[-(Upre[1,0])/(((Upre[0,0])**2.+(Upre[1,0])**2.)**(3./2.))]])
    U = Upre + DeltaT*Fpre
    t=t+DeltaT
    Orbita2D = np.append(Orbita2D,U,axis = 1)
    Upre = U
#para hacer la gráfica, desgloso las filas 1 y 2 de U como xs e ys
Xorb = Orbita2D[0,:]
Yorb = Orbita2D[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()
#
#método 2 , RK orden 4
t0=0
t=t0
#defino los vectores que voy a usar
Orbita2D = np.array([[0.],[0.],[0.],[0.]])
Orbita2D.shape = (4,1)
U = np.array([[0.],[0.],[0.],[0.]])
U.shape = (4,1)
Upre = np.array([[0.],[0.],[0.],[0.]])
Upre.shape = (4,1)
U1 = np.array([[0.],[0.],[0.],[0.]])
U1.shape = (4,1)
U2 = np.array([[0.],[0.],[0.],[0.]])
U2.shape = (4,1)
U3 = np.array([[0.],[0.],[0.],[0.]])
U3.shape = (4,1)
k1 = np.array([[0.],[0.],[0.],[0.]])
k1.shape = (4,1)
k2 = np.array([[0.],[0.],[0.],[0.]])
k2.shape = (4,1)
k3 = np.array([[0.],[0.],[0.],[0.]])
k3.shape = (4,1)
k4 = np.array([[0.],[0.],[0.],[0.]])
k4.shape = (4,1)
k0 = np.array([[0.],[0.],[0.],[0.]])
k0.shape = (4,1)
#
Upre=U0
Orbita2D=U0
while (t < tf):
    k1 = np.array([[Upre[2,0]],[Upre[3,0]],[-(Upre[0,0])/(((Upre[0,0])**2+(Upre[1,0])**2.)**(3./2.))],[-(Upre[1,0])/(((Upre[0,0])**2.+(Upre[1,0])**2.)**(3./2.))]])
    U2 = Upre + k1*DeltaT/2
    k2 = np.array([[U2[2,0]],[U2[3,0]],[-(U2[0,0])/(((U2[0,0])**2+(U2[1,0])**2.)**(3./2.))],[-(U2[1,0])/(((U2[0,0])**2.+(U2[1,0])**2.)**(3./2.))]])
    U3 = Upre + k2*DeltaT/2
    k3 = np.array([[U3[2,0]],[U3[3,0]],[-(U3[0,0])/(((U3[0,0])**2+(U3[1,0])**2.)**(3./2.))],[-(U3[1,0])/(((U3[0,0])**2.+(U3[1,0])**2.)**(3./2.))]])
    U4 = Upre + k3*DeltaT
    k4 = np.array([[U4[2,0]],[U4[3,0]],[-(U4[0,0])/(((U4[0,0])**2+(U4[1,0])**2.)**(3./2.))],[-(U4[1,0])/(((U4[0,0])**2.+(U4[1,0])**2.)**(3./2.))]])
    k0 = (1/6.)*(k1 + 2.*k2 + 2.*k3 + k4)
    U = Upre + DeltaT*k0 
    t=t+DeltaT
    Orbita2D = np.append(Orbita2D,U,axis = 1)
    Upre = U
#plot Rk4 , al compilar hay que cerrar la grafica 1 para que aparezca la nueva grafica
Xorb = Orbita2D[0,:]
Yorb = Orbita2D[1,:]
fig, ax = plt.subplots()
ax.plot(Xorb,Yorb)
plt.show()
