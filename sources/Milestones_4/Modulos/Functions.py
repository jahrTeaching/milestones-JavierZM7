from numpy import array

def Kepler(U, t): 

    x = U[0]
    y = U[1] 
    vx = U[2] 
    vy = U[3]
    d = (x**2+y**2)**1.5

    return  array( [ vx, vy, -x/d, -y/d ] )

def Oscilator_lin(U, t): 

    x = U[0]
    vx = U[1] 

    return  array( [ vx,-x ] )






