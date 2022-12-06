from Modulos.Newton_module import Newton
 
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
        return x - U - F(x,t+dt)*dt
    #return fsolve(iter, U)
    return Newton(iter, U)

def Crank_Nicolson(U,F,dt,t):
    def iter(x):
        return x - U -(F(x,t)+F(U,t))*dt/2
    return Newton(iter,U)