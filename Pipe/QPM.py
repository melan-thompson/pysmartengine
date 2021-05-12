import sys
sys.path.append("..")
def openEnd(p, p0=1.e5, T0=300, R=None, k=None):
    from GasProperty import k_Justi, Rg
    from math import sqrt
    if R is None:
        R = Rg()
    if k is None:
        k = k_Justi(T0)
    return sqrt(2 * k * R * T0 / (k - 1) * (1 - pow(p / p0, (k - 1) / k)))


def openEnd2(u, p0=1.e5, T0=300, R=None, k=None):
    from GasProperty import k_Justi, Rg
    from math import sqrt
    if R is None:
        R = Rg()
    if k is None:
        k = k_Justi(T0)
    temp = 2 * k * R * T0 / (k - 1)
    temp2 = 1 - u ** 2 / temp
    return p0 * pow(temp2, k / (k - 1))


def nozzle(p, p0=1.e5, T0=300, phi=0.5, R=None, k=None):
    from GasProperty import k_Justi, Rg
    from math import sqrt
    if R is None:
        R = Rg()
    if k is None:
        k = k_Justi(T0)
    temp = 2 * k * R * T0 / (k - 1)
    temp2 = pow(p / p0, (k - 1) / k) - 1
    temp3 = 1 / pow(phi, 2) * pow(p / p0, 2. / k) - 1

    return sqrt(temp * temp2 / temp3)
    # return sqrt(temp * temp2 * temp3)

def nozzleUn(pi,k=1.4,phi=0.5):
    temp=2/(k-1.)*(pow(pi,2)-1)
    temp2=1/pow(phi,2)*pow(pi,4./(k-1))-1
    return temp/temp2

def eqnozzle(p, u, p0=1.e5, T0=300, phi=1, R=None, k=None):
    from GasProperty import k_Justi, Rg
    from math import sqrt
    if R is None:
        R = Rg()
    if k is None:
        k = k_Justi(T0)
    temp = 2 * k * R * T0 / (k - 1)
    temp2 = pow(p / p0, (k - 1) / k) - 1
    temp3 = 1 / pow(phi, 2) * pow(p / p0, 2. / k) - 1
    # return u ** 2 * temp3 - temp * temp2
    return u ** 2 - temp * temp2*temp3


def nozzle2(u, p0=1.e5, T0=300, phi=0.9, R=None, k=None):
    from GasProperty import k_Justi, Rg
    if R is None:
        R = Rg()
    if k is None:
        k = k_Justi(T0)
    pinit = 1.1 * p0

    def fun(p):
        # return nozzle(p,p0,T0,phi,R,k)-u
        return eqnozzle(p, u, p0, T0, phi, R, k)

    h = p0 / 1.e4
    # print(fun(pinit))
    while abs(fun(pinit)) > 1.e-5:
        # print(fun(pinit))
        # print(fun(pinit+h)-fun(pinit-h))
        pinit -= fun(pinit) / ((fun(pinit + h) - fun(pinit - h)) / 2. / h)

    return pinit
