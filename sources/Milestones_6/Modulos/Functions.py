from numpy import array,size,zeros,rint
from numpy.linalg import    norm
from Modulos.CR3BPF import CR3BP
def Kepler(U, t): 

    x = U[0]
    y = U[1] 
    vx = U[2] 
    vy = U[3]
    d = (x**2+y**2)**1.5

    return  array( [ vx, vy, -x/d, -y/d ] )

def Kepler3D(U, t): 

    x = U[0]
    y = U[1] 
    z = U[2]
    vx = U[3] 
    vy = U[4]
    vz = U[5]
    d = (x**2+y**2+z**2)**1.5

    return  array( [ vx, vy, vz, -x/d, -y/d,  -z/d ] )

def Oscilator_lin(U, t): 

    x = U[0]
    vx = U[1] 

    return  array( [ vx,-x ] )


def KeplerNbody(U, t): 
    
    Nb = (size(U)/3)/2    
    Nb = int(Nb)
   
    F = zeros([2*Nb,3])
    for i in range(0,Nb): 
        F[i,:] = U[Nb+i,:]
    
    a0 = zeros([Nb,3])

    for i in range(Nb):
      for j in range(Nb):
        if j != i:
            d = U[j,:] - U[i,:] 
            a0[i,:] =  a0[i,:] + d[:]/norm(d)**3
                
    for i in range(0,Nb): 
        F[Nb+i,:] = a0[i,:]
    
    return  F

def MoonEarthCR3BP(U,t):

 mu = 1.2151e-2 #Tierra-Luna

 return CR3BP(U, mu)

def EarthSunCR3BP(U,t):

 mu = 3.0039e-7 #Sol-Tierra

 return CR3BP(U, mu)



