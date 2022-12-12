from numpy.linalg import norm, eig
from numpy import array, zeros, sqrt, size
from scipy.optimize import fsolve

def CR3BP(U, mu):
    
    p1 = sqrt((U[0]+mu)**2 + U[1]**2)
    p2 = sqrt((U[0]-1+mu)**2 + U[1]**2)

    dvdt_x = 2*U[3] + U[0] -(1-mu)*(U[0]+mu)/(p1**3) - mu*(U[0]-1+mu)/(p2**3)
    dvdt_y = -2*U[2] + U[1] -(1-mu)*U[1]/(p1**3) - mu*U[1]/(p2**3)

    return array([U[2], U[3], dvdt_x, dvdt_y])

def Lpoints(U0, Np, mu):
    LP = zeros([Np,2])

    def F(Y):
        X = zeros(4)
        X[0:2] = Y
        X[2:4] = 0
        FX = CR3BP(X, mu)
        return FX[2:4]
   
    for i in range(Np):
        LP[i,:] = fsolve(F, U0[i,0:2])

    return LP

def Jacobian(F, U):
	N = size(U)
	Jac = zeros([N,N])
	t = 1e-10

	for i in range(N):
		xj = zeros(N)
		xj[i] = t
		Jac[:,i] = (F(U + xj) - F(U - xj))/(2*t)
	return Jac

def StabLP(U0, mu):

    def F(Y):
        return CR3BP(Y, mu)

    A = Jacobian(F, U0)
    values, vectors = eig(A)

    return values