from Modulos.Scheme_module import RK4,Euler,Inverse_Euler,Crank_Nicolson, Leap_frog
from numpy import array, zeros, size, absolute, meshgrid,sqrt

def Stab_Regions(x,y,Scheme):  
    #Z = meshgrid(x,y)                         # meshgrid no permite complejos
    Z = zeros([size(x),size(x)],dtype=complex) # size(x) = size(y)
    for i in range(size(x)):                   # crea una malla XY
        for j in range(size(x)):
            #Z[size(x)-1-j,i] = complex(x[i],y[j])
             Z[j,i] = complex(x[i],y[j])

    #return absolute(array(Pol_scheme(Scheme, Z)))
    return absolute(Pol_scheme(Scheme, Z))



""" def Pol_scheme(Scheme, Z):
    if Scheme == Euler:
        r = 1 + Z
    elif Scheme == Inverse_Euler:
        r = 1/(1-Z)
    elif Scheme == Crank_Nicolson:
        r = (1+Z/2)/(1-Z/2)
    elif Scheme == RK4:
        r = 1 + Z + (Z**2)/2 + (Z**3)/(6) + (Z**4)/(24)
    elif Scheme == Leap_frog:
        r = sqrt(1) 

    return r """

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
