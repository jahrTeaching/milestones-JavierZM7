from numpy import zeros, matmul, size
from numpy.linalg import norm

def RK(order, U1, t, dt, F):
    orders, Ns, a, b, bs, c = ButcherHeunEuler()
    k = zeros([Ns, len(U1)])
    k[0,:] = F(U1, t + c[0]*dt)

    if order == 1: 
        for i in range(1,Ns):
            Up = U1
            for j in range(i):
                Up = Up + dt*a[i,j]*k[j,:]
            k[i,:] = F(Up, t + c[i]*dt)
        U2 = U1 + dt*matmul(b,k)

    elif order == 2:
        for i in range(1,Ns):
            Up = U1
            for j in range(i):
                Up = Up + dt*a[i,j]*k[j,:]
            k[i,:] = F(Up, t + c[i]*dt)
        U2 = U1 + dt*matmul(bs,k)

    return U2

def StepSize(dU, tol, dt, orders): 
    error = norm(dU)
    if error > tol:
        step_size = dt*(tol/error)**(1/(orders+1))
    else:
        step_size = dt
    return step_size

def ButcherHeunEuler():
    orders = [2,1]
    Ns = 2 

    a = zeros([Ns,Ns-1])
    b = zeros([Ns])
    bs = zeros([Ns])
    c = zeros([Ns])

    c = [0., 1.]
    a[0,:] = [0]
    a[1,:] = [  1. ]
    b[:]  = [1./2,1./2]
    bs[:] = [1.,0.]
    return orders, Ns, a, b, bs, c