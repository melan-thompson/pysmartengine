from FluidProperties.GasProperty import *
from ArrayTable import ArrayTable


class Node:
    def __init__(self):
        self.p = self.T = self.rho = self.u = None
        self.e = self.e0 = self.k = None
        self.U = None

    def Uinit(self, u, p, T, AFAK=1.e8, A=1):
        # from GasProperty import Rg
        # from GasProperty import k_Justi
        self.p = p
        self.T = T
        self.rho = p / (Rg(AFAK) * T)
        self.u = u

        self.k = 1.4  # k_Justi(self.T, AFAK)
        self.e = self.p / self.rho / (self.k - 1.)
        self.e0 = self.e + self.u ** 2 / 2.
        self.h0 = self.e0 + self.p / self.rho
        self.a = pow(self.k * self.p / self.rho, 1. / 2.)

        from numpy import array
        self.U = array([self.rho * A, self.rho * self.u * A,
                        (self.p / (self.k - 1) + 1. / 2. * self.rho * self.u ** 2) * A]).transpose()

    def characLine(self):
        self.lamda = self.a + (self.k - 1.) / 2. * self.u
        self.beta = self.a - (self.k - 1.) / 2. * self.u
        self.lamdaDirection = (self.k + 1) / 2. / (self.k - 1) * self.lamda - (3. - self.k) / 2. / (
                self.k - 1.) * self.beta
        self.betaDirection = (3. - self.k) / 2. / (self.k - 1.) * self.lamda - (self.k + 1) / 2. / (
                self.k - 1) * self.beta

    def solve(self, U=None, A=1, AFAK=1.e8):
        # from GasProperty import Rg
        # from GasProperty import k_Justi
        if U is None:
            U = self.U
        else:
            self.U = U

        self.rho = U[0] / A
        self.u = U[1] / U[0]
        self.e = U[2] / U[0] - 1. / 2. * pow(U[1] / U[0], 2)
        T = self.e * (1.4 - 1.) / Rg(AFAK)

        # def kite(T):
        #     return self.e - Rg(AFAK) * T / (k_Justi(T, AFAK) - 1)
        #
        # h = 1.
        # while abs(kite(T)) > 1.e-5:
        #     T -= kite(T) / ((kite(T + h) - kite(T - h)) / 2. / h)
        self.T = T
        self.k = k_Justi(T, AFAK)
        self.p = Rg(AFAK) * T * self.rho
        # print("T={}".format(self.T))
        # print("p={}".format(self.p))
        # print("r={}".format(self.rho))

    def F(self, U=None):
        if U is None:
            U = self.U
        from numpy import array
        temp = pow(U[1], 2) / U[0]
        self.F = array([U[1], temp + (self.k - 1) * (U[2] - temp / 2),
                        U[1] * (U[2] + (self.k - 1) * (U[2] - temp / 2)) / U[0]]).transpose()
        return self.F

    def Jaccobi(self, U=None):
        if U is None:
            self.solve()
            U = self.U
            k = self.k
        else:
            self.U = U
            self.solve()
            k = self.k
        import numpy as np
        result = np.zeros([3, 3])
        result[0][1] = 1.
        result[1][0] = (k - 3.) / 2. * (U[1] ** 2) / (U[0] ** 2)
        result[1][1] = (3. - k) * U[1] / U[0]
        result[1][2] = k - 1.
        result[2][0] = -k * U[2] * U[1] / (U[0] ** 2) + (k - 1.) * (U[1] ** 3) / (U[0] ** 3)
        result[2][1] = k * U[2] / U[0] - (k - 1.) * 3. / 2. * U[1] ** 2 / U[0] ** 2
        result[2][2] = k * U[1] / U[0]
        import numpy as np
        self.J = result
        self.Jeigenvalue, self.Jeigenvector = np.linalg.eig(self.J)
        return result

    def FVS(self):
        import numpy as np
        absA = np.dot(np.dot(self.Jeigenvector, np.diag(abs(self.Jeigenvalue))), np.linalg.inv(self.Jeigenvector))
        self.Jpositive = (self.J + absA) / 2.
        self.Jnagative = (self.J - absA) / 2.
        self.Fpositive = np.dot(self.Jpositive, self.U)
        self.Fnagative = np.dot(self.Jnagative, self.U)


def StegerWarmingFVS(t, left=[], right=[], xlim=None, numberOfNode=200, A=1):
    import numpy as np
    init = ArrayTable(2, 0)

    solution = AnalyticalSolution(t, left, right)
    if xlim is None:
        mindata = [min(solution.table[i].data) for i in range(solution.col)]
        maxdata = [max(solution.table[i].data) for i in range(solution.col)]
        xlim = [1.5 * mindata[0], 1.5 * maxdata[0]]
        print(xlim)
    deltax = (xlim[1] - xlim[0]) / numberOfNode

    for i in np.arange(xlim[0], xlim[1], deltax):
        Nodeex = Node()
        if i < 0:
            Nodeex.Uinit(left[0], left[2], left[1])
        else:
            Nodeex.Uinit(right[0], right[2], right[1])
        Nodeex.Jaccobi()
        Nodeex.FVS()
        init.append([i, Nodeex])

    thisstep = init
    tnow = 0
    Corant = 0.9

    while tnow < t:
        laststep = thisstep
        mm = 0
        for i in range(laststep.row):
            if max(abs(laststep.table[1].data[i].Jeigenvalue)) > mm:
                mm = max(abs(laststep.table[1].data[i].Jeigenvalue))

        deltat = deltax * Corant / mm
        print(deltat)
        thisstep = ArrayTable(2, 0)
        thisstep.append([laststep.table[0].data[0], laststep.table[1].data[0]])
        for i in range(1, laststep.row - 1):
            Uthis = laststep.table[1].data[i].U - deltat / deltax * (
                    (laststep.table[1].data[i + 1].Fnagative - laststep.table[1].data[i].Fnagative) + (
                    laststep.table[1].data[i].Fpositive - laststep.table[1].data[i - 1].Fpositive))
            Nodethis = Node()
            Nodethis.solve(Uthis)
            Nodethis.Jaccobi()
            Nodethis.FVS()
            thisstep.append([laststep.table[0].data[i], Nodethis])
        thisstep.append([laststep.table[0].data[-1], laststep.table[1].data[-1]])

        tnow += deltat
        print("tnow={}".format(tnow))
        print("Row of table {}".format(thisstep.row))

    result = ArrayTable(5, 0)
    result.setTableHeader(["x", "Velocity", "Density", "Pressure", "Temperature"])
    result.setTableUnit(["m", "m/s", "kg/m^3", "Pa", "K"])
    for i in range(thisstep.row):
        result.append([thisstep.table[0].data[i], thisstep.table[1].data[i].u, thisstep.table[1].data[i].rho,
                       thisstep.table[1].data[i].p,
                       thisstep.table[1].data[i].T])
    return result


def laxWendroff1step(t, left=[], right=[], xlim=None, numberOfNode=500, A=1):
    import numpy as np
    init = ArrayTable(2, 0)

    solution = AnalyticalSolution(t, left, right)
    if xlim is None:
        mindata = [min(solution.table[i].data) for i in range(solution.col)]
        maxdata = [max(solution.table[i].data) for i in range(solution.col)]
        xlim = [1.5 * mindata[0], 1.5 * maxdata[0]]
        print(xlim)
    deltax = (xlim[1] - xlim[0]) / numberOfNode

    for i in np.arange(xlim[0], xlim[1], deltax):
        Nodeex = Node()
        if i < 0:
            Nodeex.Uinit(left[0], left[2], left[1])
        else:
            Nodeex.Uinit(right[0], right[2], right[1])
        Nodeex.Jaccobi()
        Nodeex.F()
        init.append([i, Nodeex])

    thisstep = init
    tnow = 0
    Corant = 0.9

    while tnow < t:
        laststep = thisstep
        mm = 0
        for i in range(laststep.row):
            if max(abs(laststep.table[1].data[i].Jeigenvalue)) > mm:
                mm = max(abs(laststep.table[1].data[i].Jeigenvalue))

        deltat = deltax * Corant / mm
        print(deltat)
        thisstep = ArrayTable(2, 0)
        thisstep.append([laststep.table[0].data[0], laststep.table[1].data[0]])
        for i in range(1, laststep.row - 1):
            Aright=1./2.*(laststep.table[1].data[i].J+laststep.table[1].data[i+1].J)
            Aleft=1./2.*(laststep.table[1].data[i].J+laststep.table[1].data[i-1].J)
            Uthis = laststep.table[1].data[i].U\
                    - deltat /2./ deltax * (laststep.table[1].data[i + 1].F - laststep.table[1].data[i-1].F)\
                    + deltat**2 /2./ deltax**2*(np.dot(Aright,(laststep.table[1].data[i+1].F - laststep.table[1].data[i].F))
                                                -np.dot(Aleft,(laststep.table[1].data[i].F - laststep.table[1].data[i-1].F)))

            Nodethis = Node()
            Nodethis.solve(Uthis)
            Nodethis.Jaccobi()
            Nodethis.F()
            thisstep.append([laststep.table[0].data[i], Nodethis])
        thisstep.append([laststep.table[0].data[-1], laststep.table[1].data[-1]])

        tnow += deltat
        print("tnow={}".format(tnow))
        print("Row of table {}".format(thisstep.row))

    result = ArrayTable(5, 0)
    result.setTableHeader(["x", "Velocity", "Density", "Pressure", "Temperature"])
    result.setTableUnit(["m", "m/s", "kg/m^3", "Pa", "K"])
    for i in range(thisstep.row):
        result.append([thisstep.table[0].data[i], thisstep.table[1].data[i].u, thisstep.table[1].data[i].rho,
                       thisstep.table[1].data[i].p,
                       thisstep.table[1].data[i].T])
    return result



# def FVS(A):
#     import numpy as np
#     eva, evec = np.linalg.eig(A)
#     absA = np.dot(np.dot(evec, np.diag(abs(eva))), np.linalg.inv(evec))
#     Apositive = (A + absA) / 2
#     Anagative = (A - absA) / 2
#     return Apositive, Anagative


def AnalyticalSolution(t, left=[], right=[], xlim=None):
    import math
    u1 = left[0]
    u2 = right[0]
    # r1 = left[1];
    # r2 = right[1]
    T1 = left[1]
    T2 = right[1]
    p1 = left[2]
    p2 = right[2]
    # T1 = p1 / r1 / Rg();
    # T2 = p2 / r2 / Rg()
    r1 = p1 / Rg() / T1
    r2 = p2 / Rg() / T2
    k1 = k_Justi(T1)
    k2 = k_Justi(T2)
    a1 = math.sqrt(k1 * p1 / r1)
    a2 = math.sqrt(k2 * p2 / r2)

    # print("-" * 100)
    # print("{0:^50}|{1:^50}".format("ul=%.5gm/s" % u1, "ur=%.5gm/s" % u2))
    # print("{0:^50}|{1:^50}".format("rhol=%.5gkg/m^3" % r1, "rhor=%.5gkg/m^3" % r2))
    # print("{0:^50}|{1:^50}".format("pl=%sPa" % p1, "pr=%sPa" % p2))
    # print("{0:^50}|{1:^50}".format("Tl=%sK" % T1, "Tr=%sK" % T2))
    # print("-" * 100)

    def function(pbar, pj, rhoj, gamma=1.4):
        if pbar == pj:
            return 0
        elif pbar > pj:
            import math
            aj = math.sqrt(gamma * pj / rhoj)
            return (pbar - pj) / rhoj / aj / math.sqrt((gamma + 1) / 2. / gamma * pbar / pj + (gamma - 1) / 2. / gamma)
        elif pbar < pj:
            import math
            aj = math.sqrt(gamma * pj / rhoj)
            return 2. * aj / (gamma - 1.) * (pow(pbar / pj, (gamma - 1) / 2. / gamma) - 1)

    def F(p):
        return function(p, p1, r1, k1) + function(p, p2, r2, k2)

    pbar = (p1 + p2) / 2
    h = 10.
    while abs(u1 - u2 - F(pbar)) > 1.e-5:
        pbar -= (u1 - u2 - F(pbar)) / ((-F(pbar + h) + F(pbar - h)) / 2. / h)
    # print(pbar)
    ubar = (u1 + u2 - function(pbar, p1, r1, k1) + function(pbar, p2, r2, k2)) / 2.
    # print(ubar)

    A1 = r1 * a1 * math.sqrt((k1 + 1.) / (2. * k1) * pbar / p1 + (k1 - 1) / 2. / k1)
    A2 = r2 * a2 * math.sqrt((k2 + 1.) / (2. * k2) * pbar / p2 + (k2 - 1) / 2. / k2)

    if pbar >= p1:
        Z1 = u1 - A1 / r1
    else:
        Z1 = u1 - a1

    if pbar >= p2:
        Z2 = u2 + A2 / r2
    else:
        Z2 = u2 + a2

    if pbar >= p1:
        Z1star = Z1
    elif pbar > 0:
        a1star = a1 + (k1 - 1) / 2. * (u1 - ubar)
        Z1star = ubar - a1star
    else:
        Z1star = u1 - 2. * a1 / (k1 - 1)

    if pbar >= p2:
        Z2star = Z2
    elif pbar > 0:
        a2star = a2 - (k2 - 1) / 2. * (u2 - ubar)
        Z2star = ubar + a2star
    else:
        Z2star = u2 + 2 * a2 / (k2 - 1)

    r1bar = r1 * A1 / (A1 - r1 * (u1 - ubar))
    r2bar = r2 * A2 / (A2 + r2 * (u2 - ubar))

    # print("Z1=", Z1)
    # print("Z1star=", Z1star)
    # print("ubar=", ubar)
    # print("Z2star=", Z2star)
    # print("Z2=", Z2)

    if xlim is None:
        Zlaggest = abs(Z1) if abs(Z1) > abs(Z2) else abs(Z2)
        xlim = [-1.5 * t * Zlaggest, 1.5 * t * Zlaggest]

    # print("Will draw air state between [{},{}]".format(xlim[0],xlim[1]))

    def fun(t, x):
        if x / t > Z2:
            return u2, r2, p2
        elif x / t < Z1:
            return u1, r1, p1
        elif Z1star < x / t < ubar:
            return ubar, r1bar, pbar
        elif ubar < x / t < Z2star:
            return ubar, r2bar, pbar
        elif Z1 < x / t < Z1star:
            a = (k1 - 1) / (k1 + 1) * (u1 - x / t) + 2 * a1 / (k1 + 1)
            u = x / t + a
            p = p1 * pow(a / a1, 2 * k1 / (k1 - 1))
            r = k1 * p / a / a
            return u, r, p
        elif Z2star < x / t < Z2:
            a = (k2 - 1) / (k2 + 1) * (x / t - u2) + 2 * a2 / (k2 - 1)
            u = x / t - a
            p = p2 * pow(a / a2, 2 * k2 / (k2 - 1))
            r = k2 * p / a / a
            return u, r, p

    result = ArrayTable(5, 0)
    result.setTableHeader(["x", "Velocity", "Density", "Pressure", "Temperature"])
    result.setTableUnit(["m", "m/s", "kg/m^3", "Pa", "K"])
    import numpy as np
    step = (xlim[1] - xlim[0]) / 1000.
    xx = np.arange(xlim[0], xlim[1], step)
    for each in xx:
        u, r, p = fun(t, each)
        result.append([each, u, r, p, p / r / Rg()])
    return result


class ShockTube:
    def __init__(self, left, right):
        self.left = left
        self.right = right
        import math
        u1 = left[0]
        u2 = right[0]
        # r1 = left[1];
        # r2 = right[1]
        T1 = left[1]
        T2 = right[1]
        p1 = left[2]
        p2 = right[2]
        # T1 = p1 / r1 / Rg();
        # T2 = p2 / r2 / Rg()
        r1 = p1 / Rg() / T1
        r2 = p2 / Rg() / T2
        k1 = k_Justi(T1)
        k2 = k_Justi(T2)
        a1 = math.sqrt(k1 * p1 / r1)
        a2 = math.sqrt(k2 * p2 / r2)

        print("-" * 100)
        print("{0:^50}|{1:^50}".format("ul=%.5gm/s" % u1, "ur=%.5gm/s" % u2))
        print("{0:^50}|{1:^50}".format("rhol=%.5gkg/m^3" % r1, "rhor=%.5gkg/m^3" % r2))
        print("{0:^50}|{1:^50}".format("pl=%sPa" % p1, "pr=%sPa" % p2))
        print("{0:^50}|{1:^50}".format("Tl=%sK" % T1, "Tr=%sK" % T2))
        print("-" * 100)

    def AnalyticalSolutionAnimation(self, t):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        result = AnalyticalSolution(t, self.left, self.right)
        mindata = [min(result.table[i].data) for i in range(result.col)]
        maxdata = [max(result.table[i].data) for i in range(result.col)]
        fig = plt.figure(0, figsize=(10, 10))
        ax = list()
        ln = list()
        for i in range(4):
            ax.append(fig.add_subplot(4, 1, i + 1))
            plt.ylabel(result.table[i + 1].ColName + "(" + result.table[i + 1].ColUnit + ")")
            temp, = ax[i].plot([], [], 'b-', animated=False)
            ln.append(temp)
        plt.xlabel("x(m)")
        plt.tight_layout()

        def init():
            for i in range(4):
                ax[i].set_xlim(mindata[0], maxdata[0])
                ax[i].set_ylim(0.8 * mindata[i + 1], 1.05 * maxdata[i + 1])

            return ln[0], ln[1], ln[2], ln[3]

        def update(t):
            solution = AnalyticalSolution(t, self.left, self.right, [mindata[0], maxdata[0]])
            for i in range(4):
                ln[i].set_data(solution.table[0].data, solution.table[i + 1].data)

            return ln[0], ln[1], ln[2], ln[3]

        ani = FuncAnimation(fig, update, frames=np.linspace(1.e-12, t, 90), interval=100,
                            init_func=init, blit=True, repeat=False)
        ani.save("shocktube.gif", writer="pillow")
        plt.show()

    # if pbar<0:
    #     print("Left and right waves will all be rarefaction waves and there exists vacuum between them")
    #     Z1=u1-a1
    #     Z1star=u1+2*a1/(k1-1)
    #     Z2=u2+a2
    #     Z2star=u2-2*a2/(k2-1)
    #     def fun(t,x):
    #         if x/t<Z1: #左行波波前
    #             return u1,r1,p1
    #         if Z1star>x/t>Z1: #左行稀疏波波内
    #             a=(k1-1)/(k1+1)*(u1-x/t)+2*a1/(k1+1)
    #             u=x/t+a
    #             p=p1*pow(a/a1,2*k1/(k1-1))
    #             r=k1*p/a/a
    #             return u,r,p
    #         if x/t>Z2:  #右行波波前
    #             return u2,r2,p2
    #         if Z2>x/t>Z2star: #右行稀疏波波内
    #             a=(k2-1)/(k2+1)*(x/t-u2)+2*a2/(k2-1)
    #             u=x/t-a
    #             p=p2*pow(a/a2,2*k2/(k2-1))
    #             r=k2*p/a/a
    #             return u,r,p

    # if p2 >= p1:
    #     if (u1 - u2) >= F(p2):
    #         print("Left and right waves will all be shock waves")
    #         a1 = math.sqrt(k1 * r1 / p1)
    #         A1 = r1 * a1 * math.sqrt((k1 + 1) / 2. / k1 * pbar / p1 + (k1 - 1) / 2. / k1)
    #         Z1 = u1 - A1 / r1
    #         rho1bar = r1 * A1 / (A1 - r1 * (u1 - ubar))
    #         a2=math.sqrt(k2*r2/p2)
    #         A2=r2*a2*math.sqrt((k2+1)/2./k2*pbar/p2+(k2-1)/2./k2)
    #         Z2=u2+A2/r2
    #         rho2bar=r2*A2/(A2+r2*(u2-ubar))
    #         def fun(t,x):
    #             if Z1>x/t:
    #                 return u1,r1,p1
    #             elif Z1 < x/t < ubar:
    #                 return ubar,rho1bar,pbar
    #             elif x/t>Z2:
    #                 return u2,r2,p2
    #             elif Z2 > x/t > ubar:
    #                 return ubar,rho2bar,pbar
    #     elif (u1 - u2) >= F(p1):
    #         print("left wave will be shock wave and right wave will be rarefaction wave")
    #
    #     elif (u1 - u2) >= F(0.):
    #         print("Left and right waves will all be rarefaction waves")
    #     else:
    #         print("Left and right waves will all be rarefaction waves and vacuum will between them")
    # else:
    #     if (u1 - u2) >= F(p1):
    #         print("Left and right waves will all be shock waves")
    #     elif (u1 - u2) >= F(p2):
    #         print("right wave will be shock wave and left wave will be rarefaction wave")
    #     elif (u1 - u2) >= F(0.):
    #         print("Left and right waves will all be rarefaction waves")
    #     else:
    #         print("Left and right waves will all be rarefaction waves and vacuum will between them")
    def StegerWarmingFVSAnimation(self):
        pass


if __name__=="__main__":
    ShockTube([1,1200,0.5e6],[0,300,0.1e6]).AnalyticalSolutionAnimation(0.6e-3)



class pipe1D:
    def __init__(self, pu, pd, Tu, Td):
        self.pu = pu
        self.pd = pd

    def boundaryCondition(self):
        pass
