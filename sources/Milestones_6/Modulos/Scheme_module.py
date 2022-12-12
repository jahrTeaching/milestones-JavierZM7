from scipy.optimize import fsolve
from Modulos.RKEMB import ButcherHeunEuler,RK,StepSize

def Euler(U,F,dt,t):

    return U + dt*F(U,t)

def RK4(U,F,dt,t):
    k1=F(U,t)
    k2=F(U + dt*k1/2 , t)
    k3=F(U + dt*k2/2 , t)
    k4=F(U + dt*k3 , t)
    k=(1/6)*(k1 + 2*k2 + 2*k3 + k4)

    return U + dt*k

def Inverse_Euler(U,F,dt,t):

    def iter(x):
        return x - U - F(x,t)*dt
    #return fsolve(iter, U)
    return fsolve(iter, U)

def Crank_Nicolson(U,F,dt,t):
    def iter(x):
        return x - U -(F(x,t)+F(U,t))*dt/2
    return fsolve(iter,U)

def Leap_frog(U1,U2,F,dt,t):
    
    return U1 + 2*dt*F(U2,t)

def RKEmb(U, F, dt, t):

    tol = 1e-9
    orders, Ns, a, b, bs, c = ButcherHeunEuler()
    est1 = RK(1, U, t, dt, F) 
    est2 = RK(2, U, t, dt, F) 
    h = min(dt, StepSize(est1-est2, tol, dt,  min(orders)))
    N_n = int(dt/h)+1
    n_dt = dt/N_n
    est1 = U
    est2 = U

    for i in range(N_n):
        time = t + i*dt/int(N_n)
        est1 = est2
        est2 = RK(1, est1, time, n_dt, F)

    return est2