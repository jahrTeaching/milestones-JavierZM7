from numpy import linspace,size,zeros,array,linalg,shape,log10
from Modulos.Cauchy_module import Cauchy
from Modulos.Scheme_module import Crank_Nicolson,RK4 

def Error_R(Function,U0,linspace_t,scheme):
    N = size(linspace_t)
    linspace_t2 = linspace(0,linspace_t[N-1],2*N-1)

    if scheme == RK4:
        q = 4
    elif scheme == Crank_Nicolson:
        q = 2
    else:
        q = 1

    E = zeros([size(U0),N])
    U2n = Cauchy(Function, U0, linspace_t2, scheme)
    U1n = Cauchy(Function, U0, linspace_t, scheme)

    for i in range (0,N):
        E[:,i]=((U2n[:,2*i]-U1n[:,i])/(1-1/(2**q)))
        
    return E

def Convergence_R(Function,U0,linspace_t,scheme): 
    
    p = 6 #6                                                             #number of points
    Nt2 = size(linspace_t)
    U1= Cauchy(Function, U0, linspace_t, scheme)
    
    log_E = zeros(p)
    log_N = zeros(p)
    log_Eq = zeros(p)

    U2 = zeros( shape(U1) )

    for i in range(p):
        print(i)
        Nt2 = 2*Nt2
        linspace_t2 = linspace(0, linspace_t[-1],Nt2)                                    #t[-1] return the last value of the vector or list
        U2 = Cauchy(Function, U0, linspace_t2, scheme)
        log_E[i] = log10(linalg.norm(  U2[:, -1]- U1[:,-1] ) )
        log_N[i] = log10(Nt2)
        linspace_t = linspace_t2
        U1 = U2
            
    q = -((log_E[3]-log_E[2])/(log_N[3]-log_N[2]))
    cq = log10(1-(1/(2**q)))                                            #constante logaritmica de q

    cq_vector = array([cq,cq,cq,cq,cq,cq,cq,cq,cq,cq,cq,cq])                  
    for i in range(p):
        
        log_Eq[i] = log_E[i] - cq_vector[i]
            
    return log_Eq, log_N ,q

