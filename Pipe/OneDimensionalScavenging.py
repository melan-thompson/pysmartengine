from FluidProperties.GasProperty import *


def DavisFluxLimiter(r: float) -> float:
    if r <= 0:
        return 0
    else:
        return min(2 * r, 1)


def C(v):
    if v > 0.5:
        return 0.25
    elif v <= 0.5:
        return v * (1 - v)


def fluxlimeter(v, r):
    return 1. / 2. * C(v) * (1 - DavisFluxLimiter(r))


def fraction(v, D, f=0.005):
    # G = \frac{{\tau \pi Ddx}}{{\rho Adx}} = \frac{{f\left( {\frac{1}{2}\rho {u^2}} \right)\pi \sqrt {\frac{{4A}}{\pi }} }}{{\rho A}} = f{u^2}\sqrt {\frac{\pi }{A}}  = \frac{{f\dot m\left| {\dot m} \right|}}{{{A^2}{\rho ^2}}}\sqrt {\frac{\pi }{A}}
    # G = \frac{{\pi D\tau dx}}{{dm}}
    """
    单位质量流体分配到沿x轴方向的比摩擦力，
    :param v:流体速度m/s
    :param D:管道直径,m
    :param f:管道的摩擦系数，定义为f{\rm{ = }}\frac{{{\tau _0}}}{{\frac{1}{2}\rho {c^2}}}，一般取0.004-0.01
    :return:管道摩擦力
    """
    if v == 0:
        return 0
    else:
        return 2 * f * (v / abs(v)) * pow(v, 2) / D


def massFraction2Alphak(Y, L0=14.3):
    if Y > 1 or Y < 0:
        raise Exception("mass fraction Y={} out of range 0 to 1".format(Y))
    elif Y == 1:
        return 1.e8
    else:
        return 1 + (1 + L0) / L0 * Y / (1 - Y)


class Node:
    def __init__(self, r=1):
        self.p = self.T = self.rho = self.u = self.a = None
        self.e = self.e0 = self.k = self.h0 = None
        self.Y = self.w = self.alphak = None
        self.U = None
        self.r = r

    def Uinit(self, u, T, p, w, Y):
        self.alphak = massFraction2Alphak(Y)
        self.p = p
        self.T = T
        self.rho = p / Rg(self.alphak) / T
        self.u = u
        self.w = w
        self.Y = Y

        self.k = k_Justi(T, self.alphak)
        self.e = self.p / self.rho / (self.k - 1)
        self.e0 = self.e + self.u ** 2 / 2.
        self.h0 = self.e0 + self.p / self.rho
        self.a = pow(self.k * self.p / self.rho, 1. / 2.)
        # print(self.rho)
        # print(self.u)
        # print(self.rho * w * self.r)

        from numpy import array
        self.U = array(
            [self.rho, self.rho * self.u, self.rho * self.e0, self.rho * w * self.r, self.rho * self.Y]).transpose()

    def solve(self, U=None):
        if U is None:
            U = self.U
        else:
            self.U = U

        self.rho = U[0]
        self.u = U[1] / U[0]
        self.Y = U[4] / U[0]
        if self.Y > 1: self.Y = 1
        if self.Y < 0: self.Y = 0

        self.e = U[2] / U[0] - 1. / 2. * pow(self.u, 2)

        self.alphak = massFraction2Alphak(self.Y)

        self.T = self.e * (1.35 - 1) / Rg(self.alphak)

        self.k = k_Justi(self.T, self.alphak)
        self.p = (self.k - 1) * (U[2] - 1. / 2. * pow(U[1], 2) / U[0])
        self.w = U[3] / self.rho / self.r

        self.e0 = self.e + 1. / 2. * self.u ** 2
        self.h0 = self.e0 + self.p / self.rho

        self.a = pow(self.k * self.p / self.rho, 1. / 2.)

    def Flux(self, U=None):
        if U is None:
            U = self.U
        from numpy import array

        temp = pow(U[1], 2) / U[0]  # rho*u^2
        self.F = array([U[1], (3 - self.k) / 2 * temp + (self.k - 1) * U[2],
                        U[1] / U[0] * (self.k * U[2] - (self.k - 1) / 2 * temp), U[3] * U[1] / U[0],
                        U[4] * U[1] / U[0]])
        return self.F

    def Source(self, U=None):
        if U is None:
            U = self.U
        from numpy import array
        w = U[3] / U[0] / self.r
        self.S = array([0, U[0] * fraction(U[1] / U[0], self.r * 2), 0, U[0] * w ** 2 * self.r ** 2 * 0.005, 0])

    def Jaccobi(self, U=None):
        # 首先将需要用到的量求出来
        if U is None:
            self.solve()
            U = self.U
        else:
            self.U = U
            self.solve()

        k = self.k

        import numpy as np
        result = np.zeros([5, 5])
        result[0][1] = 1

        result[1][0] = (k - 3) / 2. * self.u ** 2
        result[1][1] = (3 - k) * self.u
        result[1][2] = k - 1

        result[2][0] = -k * self.u * self.e0 + (k - 1) * self.u ** 3
        result[2][1] = k * self.e0 + 3. / 2. * (1 - k) * self.u ** 2
        result[2][2] = k * self.u

        result[3][0] = -self.u * self.w * self.r
        result[3][1] = self.w * self.r
        result[3][3] = self.u

        result[4][0] = -self.u * self.Y
        result[4][1] = self.Y
        result[4][4] = self.u

        self.J = result
        self.Jeigenvalue, self.Jeigenvector = np.linalg.eig(self.J)
        return result

    def printProperties(self):
        print(
            "rho={}kg/m^3,p={}Pa,T={}K,u={}m/s,Y={},w={},a={}".format(self.rho, self.p, self.T, self.u, self.Y, self.w,
                                                                      self.a))

        print("k={},alphak={}".format(self.k, self.alphak))

        print("e={}J/kg,e0={}J/kg,h0={}J/kg".format(self.e, self.e0, self.h0))


def twoStepLaxWendroff(t, left=[], right=[], xlim=None, number_of_nodes=500, D=1):
    from ArrayTable import ArrayTable
    from ShockTube import AnalyticalSolution
    solution = AnalyticalSolution(t, left=left[0:3], right=right[0:3])
    # solution.plot(1)
    # solution.plot(2)
    # solution.plot(3)
    # solution.plot(4)

    import numpy as np
    init = ArrayTable(2, 0)
    if xlim is None:
        mindata = [min(solution.table[i].data) for i in range(solution.col)]
        maxdata = [max(solution.table[i].data) for i in range(solution.col)]
        xlim = [1.5 * mindata[0], 1.5 * maxdata[0]]
        print(xlim)
    deltax = (xlim[1] - xlim[0]) / number_of_nodes

    for i in np.arange(xlim[0], xlim[1], deltax):
        Nodeex = Node()
        if i < 0:
            Nodeex.Uinit(left[0], left[1], left[2], left[3], left[4])
        else:
            Nodeex.Uinit(right[0], right[1], right[2], right[3], right[4])
        Nodeex.Jaccobi()
        Nodeex.Flux()
        Nodeex.Source()
        init.append([i, Nodeex])

    thisstep = init
    tnow = 0
    Corant = 0.6

    while tnow < t:
        laststep = thisstep
        maxeigenvalue = 0
        for i in range(laststep.row):
            if max(abs(laststep.table[1].data[i].Jeigenvalue)) > maxeigenvalue:
                maxeigenvalue = max(abs(laststep.table[1].data[i].Jeigenvalue))

        deltat = deltax * Corant / maxeigenvalue
        print("t step={}".format(deltat))
        thisstep = ArrayTable(2, 0)

        thisstep.append([laststep.table[0].data[0], laststep.table[1].data[0]])
        thisstep.append([laststep.table[0].data[1], laststep.table[1].data[1]])

        for i in range(2, laststep.row - 2):
            Uhalfright = 1. / 2. * (
                    laststep.table[1].data[i].U + laststep.table[1].data[i + 1].U) - deltat / 2. / deltax * (
                                 laststep.table[1].data[i + 1].F - laststep.table[1].data[i].F) - deltat / 4. * (
                                 laststep.table[1].data[i + 1].S + laststep.table[1].data[i].S)
            # print(Uhalfright)
            Nodehalfright = Node()
            Nodehalfright.solve(Uhalfright)
            Nodehalfright.Flux()
            Nodehalfright.Source()

            Uhalfleft = 1. / 2. * (
                    laststep.table[1].data[i].U + laststep.table[1].data[i - 1].U) - deltat / 2. / deltax * (
                                laststep.table[1].data[i].F - laststep.table[1].data[i - 1].F) - deltat / 4. * (
                                laststep.table[1].data[i - 1].S + laststep.table[1].data[i].S)
            Nodehalfleft = Node()
            Nodehalfleft.solve(Uhalfleft)
            Nodehalfleft.Flux()
            Nodehalfleft.Source()

            # 未加人工粘性
            # Uthis = laststep.table[1].data[i].U - deltat / deltax * (Nodehalfright.F - Nodehalfleft.F) - deltat / 2. * (
            #         Nodehalfright.S + Nodehalfleft.S)

            # 加入人工耗散项
            Up_1_2 = laststep.table[1].data[i + 1].U - laststep.table[1].data[i].U
            Un_1_2 = laststep.table[1].data[i].U - laststep.table[1].data[i - 1].U

            def r(w1, w2, w3, w4=None):
                if w4 is None:
                    w4 = w3
                if np.all(w3 == 0): return 0
                return np.inner(w1, w2) / np.inner(w3, w4)

            rp_i_1 = r(laststep.table[1].data[i - 1].U - laststep.table[1].data[i - 2].U,
                       laststep.table[1].data[i].U - laststep.table[1].data[i - 1].U,
                       Un_1_2, Un_1_2)
            # print(laststep.table[1].data[i - 1].U - laststep.table[1].data[i - 2].U)
            # print(rp_i_1)
            rn_i = r(laststep.table[1].data[i].U - laststep.table[1].data[i - 1].U,
                     laststep.table[1].data[i + 1].U - laststep.table[1].data[i].U,
                     Un_1_2, Un_1_2)

            rp_i = r(laststep.table[1].data[i].U - laststep.table[1].data[i - 1].U,
                     laststep.table[1].data[i + 1].U - laststep.table[1].data[i].U,
                     Up_1_2, Up_1_2)

            rn_i_1 = r(laststep.table[1].data[i + 1].U - laststep.table[1].data[i].U,
                       laststep.table[1].data[i + 2].U - laststep.table[1].data[i + 1].U,
                       Up_1_2, Up_1_2)

            v = abs(maxeigenvalue) * deltat / deltax

            Uthis = laststep.table[1].data[i].U - deltat / deltax * (Nodehalfright.F - Nodehalfleft.F) - deltat / 2. * (
                    Nodehalfright.S + Nodehalfleft.S) + (fluxlimeter(v, rp_i) + fluxlimeter(v, rn_i_1)) * Up_1_2 - (
                            fluxlimeter(v, rp_i_1) + fluxlimeter(v, rn_i)) * Un_1_2
            Nodethis = Node()
            Nodethis.solve(Uthis)
            Nodethis.Jaccobi()
            Nodethis.Flux()
            Nodethis.Source()
            thisstep.append([laststep.table[0].data[i], Nodethis])
        thisstep.append([laststep.table[0].data[-1], laststep.table[1].data[-1]])
        thisstep.append([laststep.table[0].data[-2], laststep.table[1].data[-2]])

        tnow += deltat
        print("tnow={}".format(tnow))
        # print("Row of table {}".format(thisstep.row))

        result = ArrayTable(7, 0)
        result.setTableHeader(["x", "Velocity", "Density", "Pressure", "Temperature", "$omega$", "air mass fraction"])
        result.setTableUnit(["m", "m/s", "kg/m^3", "Pa", "K", "deg/s", "/"])
        for i in range(thisstep.row):
            result.append([thisstep.table[0].data[i], thisstep.table[1].data[i].u, thisstep.table[1].data[i].rho,
                           thisstep.table[1].data[i].p,
                           thisstep.table[1].data[i].T, thisstep.table[1].data[i].w, thisstep.table[1].data[i].Y])
        yield result


class OneDimensionalScavenge:
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def animation(self, t):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        from ShockTube import AnalyticalSolution
        # result = AnalyticalSolution(t, self.left[0:3], self.right[0:3])
        # mindata = [min(result.table[i].data) for i in range(result.col)]
        # maxdata = [max(result.table[i].data) for i in range(result.col)]
        fig = plt.figure(0, figsize=(10, 10))

        s = twoStepLaxWendroff(t, self.left, self.right)
        result = next(s)
        result = next(s)
        result = next(s)
        result = next(s)
        mindata = [min(result.table[i].data) for i in range(result.col)]
        maxdata = [max(result.table[i].data) for i in range(result.col)]
        # result.plot()

        ax = list()
        ln = list()
        for i in range(6):
            ax.append(fig.add_subplot(6, 1, i + 1))
            plt.ylabel(result.table[i + 1].ColName + "(" + result.table[i + 1].ColUnit + ")")
            temp, = ax[i].plot([], [], 'b-', animated=False)
            ln.append(temp)
        plt.xlabel("x(m)")
        plt.tight_layout()

        def init():
            for i in range(6):
                ax[i].set_xlim(0.5 * mindata[0], 0.5 * maxdata[0])
                ax[i].set_ylim(0.8 * mindata[i + 1], 1.2 * maxdata[i + 1])

            return ln[0], ln[1], ln[2], ln[3], ln[4], ln[5]

        def update(t):
            try:
                solution = next(s)
            except StopIteration:
                print("iteration stops")
                return None
            # solution.plot(1)
            for i in range(6):
                ln[i].set_data(solution.table[0].data, solution.table[i + 1].data)

            return ln[0], ln[1], ln[2], ln[3], ln[4], ln[5]

        ani = FuncAnimation(fig, update, frames=np.linspace(1.e-12, t, 10000), interval=100,
                            init_func=init, blit=True, repeat=False)
        # ani.save("withfluxlimiter.gif", writer="pillow")

        plt.show()


if __name__ == "__main__":
    result = OneDimensionalScavenge([0, 300, 5.e5, 0, 1], [0, 1000, 1.5e5, 1, 0])
    result.animation(1.e-2)
    #
    s = twoStepLaxWendroff(1.e-2, [0, 300, 5.e5, 0, 1], [0, 1000, 1.5e5, 1, 0])
    s.plot(1)
    s.plot(2)
    s.plot(3)
    s.plot(4)
    s.plot(5)
    s.plot(6)
    #
    # import numpy as np
    # import matplotlib.pyplot as plt
    # from matplotlib.animation import FuncAnimation

    # result = AnalyticalSolution(t, self.left, self.right)
    # mindata = [min(result.table[i].data) for i in range(result.col)]
    # maxdata = [max(result.table[i].data) for i in range(result.col)]
    # fig = plt.figure(0, figsize=(10, 10))
    # ax = list()
    # ln = list()
    # for i in range(4):
    #     ax.append(fig.add_subplot(4, 1, i + 1))
    #     plt.ylabel(result.table[i + 1].ColName + "(" + result.table[i + 1].ColUnit + ")")
    #     temp, = ax[i].plot([], [], 'b-', animated=False)
    #     ln.append(temp)
    # plt.xlabel("x(m)")
    # plt.tight_layout()
    #
    #
    # def init():
    #     for i in range(4):
    #         ax[i].set_xlim(mindata[0], maxdata[0])
    #         ax[i].set_ylim(0.8 * mindata[i + 1], 1.05 * maxdata[i + 1])
    #
    #     return ln[0], ln[1], ln[2], ln[3]
    #
    #
    while True:
        try:
            # next(s).plot(1)

            # next(s).plot(2)
            next(s).plot(3)
            # next(s).plot(4)
            # next(s).plot(5)
            # next(s).plot(6)
        except StopIteration:
            break

    # N = Node()
    # N.Uinit(10,300,10, 1, 1)
    # N.printProperties()
    # print(N.Flux())
    # print(N.Jaccobi())
    # print(N.Jeigenvalue)
    # print(N.Jeigenvector)
