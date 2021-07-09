class OneDimensionalScavenge:
    def __init__(self,rho,u,e,Y,omega,r):
        from torch import tensor
        U=tensor([rho,rho*u,rho*omega*r,rho*e,rho*Y])