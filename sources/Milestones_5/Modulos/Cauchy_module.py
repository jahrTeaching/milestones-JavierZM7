from Modulos.Scheme_module import Leap_frog
from numpy import zeros

def Cauchy(F,U0,lin_t,scheme):                           #aqui se construye la matriz con las n columnas de soluciones usando cualquier esquema
                                                    # se supone que U0 es dato y se inserta como un array columna [N,1], la matriz U se tiene que definir con N filas y N+1 columnas, contando con que la primera de ellas va a ser el array U0
    
    N = lin_t.shape[0]-1
    U = zeros((U0.shape[0],U0.shape[1],N+1))
    U.shape = (U0.shape[0],U0.shape[1],N+1)
    U[:,:,0] = U0
    dt = lin_t[1]-lin_t[0]
    if scheme == Leap_frog:
        U[:,:,1] = U[:,:,0] + dt*F(U[:,:,0],lin_t[0])      #Euler(U[:,0],F,lin_t[1]-lin_t[0],lin_t[0])
        for i in range(1,N):
            
            U1 = U[:,:,i-1]
            U2 = U[:,:,i]
            U[:,:,i+1] = scheme(U1,U2,F,lin_t[i+1]-lin_t[i],lin_t[i])
    else:

     for i in range(N):
        U[:,:,i+1] = scheme(U[:,:,i],F,lin_t[i+1]-lin_t[i],lin_t[i]) # los input son U por columnas, Función derivada, salto de tiempo , instante

    return U

