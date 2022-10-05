
from numpy import zeros

def Cauchy(Function,U0,linspace_t,scheme):                           #aqui se construye la matriz con las n columnas de soluciones usando cualquier esquema
                                                    # se supone que U0 es dato y se inserta como un array columna [N,1], la matriz U se tiene que definir con N filas y N+1 columnas, contando con que la primera de ellas va a ser el array U0
    
    N = linspace_t.shape[0]-1
    U = zeros((U0.shape[0],N+1))
    U.shape = (U0.shape[0],N+1)
    U[:,0] = U0
    for i in range(N):
        U[:,i+1] = scheme(U[:,i],Function,linspace_t[i+1]-linspace_t[i],linspace_t[i]) # los input son U por columnas, Función derivada, salto de tiempo , instante

    return U