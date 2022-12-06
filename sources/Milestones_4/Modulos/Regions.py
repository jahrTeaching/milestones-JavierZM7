from Modulos.Scheme_module import RK4,Euler,Inverse_Euler,Crank_Nicolson, Leap_frog
from numpy import zeros, size, absolute, meshgrid,sqrt

def Stab_Regions(x,y,Scheme):  
    #Z = meshgrid(x,y) 
    p = size(x)                                                                                                           # meshgrid no permite complejos
    py = size(y)
    if p == py:
        Z = zeros([p,p],dtype=complex)                                                                                    # size(x) = size(y)
        for i in range(p):                                                                                                # crea una malla XY
            for j in range(p):
               Z[j,i] = complex(x[i],y[j])
    else:
        print('WARNING: X partitions are not equal to Y, meshgrid incompatible with stab_Regions function')

    
    return absolute(Pol_scheme(Scheme, Z))

def Pol_scheme(Scheme, Z):
    if Scheme == Euler:
        return  1 + Z
    elif Scheme == Inverse_Euler:
        return 1/(1-Z)
    elif Scheme == Crank_Nicolson:
        return (1+Z/2)/(1-Z/2)
    elif Scheme == RK4:
        return 1 + Z + (Z**2)/2 + (Z**3)/(6) + (Z**4)/(24)
    elif Scheme == Leap_frog:
        return sqrt(1)


