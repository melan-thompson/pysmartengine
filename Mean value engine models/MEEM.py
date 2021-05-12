# engine mean value model only calculates mean values of each variables during one or more cylcles
# it does not calculate crank angle resolved parameters

def xi(theta, theta0, theta1):
    if theta < theta0 or theta > theta1:
        raise Exception("Wrong error!!")
    from math import cos, pi
    return 0.5 * (1 - cos(pi * (theta - theta0) / (theta1 - theta0)))


class IdealCylcle:
    def __init__(self, CylGeo, ValveTiming=None):
        from ArrayTable import ArrayTable
        self.CompressData = ArrayTable(5, 0)
        self.ExpanseData = ArrayTable(5, 0)
        self.Rpressure = ArrayTable(2, 0)
        self.GasExchange = ArrayTable(2, 0)
        self.CompressData.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        self.CompressData.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])
        self.ExpanseData.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        self.ExpanseData.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])
        self.data = ArrayTable(2, 0)
        self.data.setTableHeader(["Crank angle", "Cylinder pressure"])
        self.data.setTableUnit(["°CA", "Pa"])

        self.CylGeo = CylGeo
        from Valve import ValveDesign
        if ValveTiming is None:
            self.ValveTiming = ValveDesign()
        else:
            self.ValveTiming = ValveTiming

    def compress(self, ivc, pim=1.e5, Tim=300, Tr=800, xr=0.0, kc=1.3):
        from GasProperty import DieselMixture, Rg
        self.ivc = ivc
        self.pim = pivc = pim
        self.Tim = Tivc = Tim * (1 - xr) + xr * Tr
        Vivc = self.CylGeo.V(ivc)
        # mivc = pivc * Vivc / Rg() / Tivc
        mivc = pivc * self.CylGeo.V(180) /  Rg() / Tivc
        print("Intake air mass {} mg".format(mivc*1.e6))
        self.mix = DieselMixture()
        self.mix.init_With_Mc_r(mivc, xr, mivc / 14.3 / 1.1)
        print("Tivc=", Tivc)

        for i in range(ivc, ivc + 720):
            V = self.CylGeo.V(i)
            T = Tivc * pow(Vivc / V, kc - 1)
            p = pivc * pow(Vivc / V, kc)
            m = p * V / self.mix.Rg_gas(0) / T
            self.CompressData.append([i, V, p, T, m])

        from Valve import mod
        for i in range(self.CompressData.row):
            self.CompressData.table[0].data[i] = mod(self.CompressData.table[0].data[i])
        self.CompressData.doQuickSort(0)

    def Burn(self, Rp, SOC=-5, alpha=1, Hu=42700e3, L0=14.3):
        from numpy import arange
        for i in arange(self.ValveTiming.IVC, SOC):
            self.data.append([i, self.CompressData.linearInterpolate(i, 2)])

        if Rp < 0 or Rp > 1: raise Exception("premixed burn fraction value error!!")
        self.Rp = Rp
        self.L0 = L0
        self.alpha = alpha
        self.SOC = SOC
        self.mix.gf = self.mix.M_air(0) / L0 / alpha

        def fun(x, T):
            return (Hu - self.mix.u(x, T)) / (alpha * L0 + x) / self.mix.cv(x, T)

        x = 0
        self.T2 = T = self.CompressData.linearInterpolate(SOC, 3)
        self.p2 = self.CompressData.linearInterpolate(SOC, 2)

        step = Rp / 50.
        while x < Rp:
            T += fun(x, T) * step
            x += step

        print("Temperature after premixed burn {}".format(T))
        Ttemp = T
        self.p3 = self.p2 * T / self.T2
        print("Pressure after premixed burn {}".format(self.p3))

        def fun2(x, T):
            return (Hu - self.mix.h(x, T)) / (self.mix.cp(x, T) * (self.alpha * self.L0 + x))

        step2 = (1 - self.Rp) / 100.

        while x < 1:
            T += fun2(x, T) * step2
            x += step2

        print("Temperature after disfusion burn {}".format(T))
        self.V3 = self.CylGeo.V(0) * T / Ttemp
        self.EOC = self.CylGeo.getFi(self.V3)
        print("End of combustion {}".format(self.EOC))
        self.T3 = T
        # self.p3=self.mix.M_total(1)*self.mix.Rg_gas(1)*self.T3/self.V3
        # print("Pressure after burned {}".format(self.p3))

        return T

    def Expense(self, ke=1.33):
        from numpy import arange
        for i in arange(self.SOC, self.SOC + 360):
            V = self.CylGeo.V(i)
            p = self.p3 * pow(self.V3 / V, ke)
            T = self.T3 * pow(self.V3 / V, ke - 1)
            self.ExpanseData.append([i, V, p, T, self.mix.M_total(1)])

        from Valve import mod
        for i in range(self.ExpanseData.row):
            self.ExpanseData.table[0].data[i] = mod(self.ExpanseData.table[0].data[i])
        self.ExpanseData.doQuickSort(0)

        from numpy import arange
        for i in arange(self.EOC, self.ValveTiming.EVO):
            self.data.append([i, self.ExpanseData.linearInterpolate(i, 2)])

    def pressureReconstruct(self):
        from Cylinder import WibeFunction
        from numpy import arange
        for i in arange(self.SOC, self.EOC):
            temp = WibeFunction(i, self.SOC, self.EOC - self.SOC, m=1)
            p = (1 - temp) * self.CompressData.linearInterpolate(i, 2) + temp * self.ExpanseData.linearInterpolate(i, 2)
            self.Rpressure.append([i, p])

        from numpy import arange
        for i in arange(self.SOC, self.EOC):
            self.data.append([i, self.Rpressure.linearInterpolate(i, 1)])

    def gasExchange(self,pem=3.e5):
        self.pem = pem
        from numpy import arange
        evc = self.ValveTiming.EVC
        int = self.ValveTiming.int = (self.ValveTiming.EVC + self.ValveTiming.IVC) / 2.
        for i in arange(self.ValveTiming.EVC, int):
            self.data.append([i, self.pim])
        for i in arange(int, self.ivc):
            temp = xi(i, int, self.ivc)
            p = self.pim * (1 - temp) + self.CompressData.linearInterpolate(i, 2) * temp
            self.data.append([i, p])

        exh = self.ValveTiming.exh = (self.ValveTiming.EVO + self.ValveTiming.IVO) / 2.
        for i in arange(self.ValveTiming.EVO, exh):
            temp = xi(i, self.ValveTiming.EVO, exh)
            p = self.ExpanseData.linearInterpolate(i, 2) * (1 - temp) + pem * temp
            self.data.append([i, p])
        for i in arange(exh, self.ValveTiming.IVO):
            self.data.append([i, pem])

        from Valve import mod
        for i in arange(self.ValveTiming.IVO, self.ValveTiming.EVC + 720):
            temp = xi(i, self.ValveTiming.IVO, self.ValveTiming.EVC + 720)
            p = pem * (1 - temp) + self.pim * temp
            self.data.append([mod(i), p])

        # from numpy import arange
        # for i in arange(self.ValveTiming.EVC, self.ValveTiming.IVC):
        #     self.data.append([i, self.GasExchange.linearInterpolate(i, 1)])

    def plot(self):
        self.data.doQuickSort(0)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 10))
        ax.plot(self.data.table[0].data, self.data.table[1].data)
        plt.xlabel(self.data.table[0].ColName + "(" + self.data.table[0].ColUnit + ")")
        plt.ylabel(self.data.table[1].ColName + "(" + self.data.table[1].ColUnit + ")")

        ax.scatter(self.ValveTiming.IVC, self.data.linearInterpolate(self.ValveTiming.IVC, 1))
        ax.annotate('IVC %.3g $^\circ$CA' % self.ValveTiming.IVC,
                    xy=(self.ValveTiming.IVC, self.data.linearInterpolate(self.ValveTiming.IVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))
        ax.scatter(self.EOC, self.data.linearInterpolate(self.EOC, 1))
        ax.annotate('EOC %.3g $^\circ$CA' % self.EOC,
                    xy=(self.EOC, self.data.linearInterpolate(self.EOC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        index = self.data.findMaxValueIndex(1)
        maxpreCA = self.data.table[0].data[index]
        ax.scatter(maxpreCA, self.data.linearInterpolate(maxpreCA, 1))
        ax.annotate('maxium pressure %.5g bar' % (self.data.linearInterpolate(maxpreCA, 1) / 1.e5),
                    xy=(maxpreCA, self.data.linearInterpolate(maxpreCA, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        # ax.scatter(self.IVO, table.linearInterpolate(self.IVO, 1))
        # ax.annotate('IVO %.3g $^\circ$CA' % self.IVO,
        #             xy=(self.IVO, table.linearInterpolate(self.IVO, 1)), xycoords='data',
        #             xytext=(0, 10), textcoords='offset points',
        #             arrowprops=dict(arrowstyle="->"))
        #
        ax.scatter(self.ValveTiming.EVO, self.data.linearInterpolate(self.ValveTiming.EVO, 1))
        ax.annotate('EVO %.3g $^\circ$CA' % self.ValveTiming.EVO,
                    xy=(self.ValveTiming.EVO, self.data.linearInterpolate(self.ValveTiming.EVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.ValveTiming.EVC, self.data.linearInterpolate(self.ValveTiming.EVC, 1))
        ax.annotate('EVC %.3g $^\circ$CA' % self.ValveTiming.EVC,
                    xy=(self.ValveTiming.EVC, self.data.linearInterpolate(self.ValveTiming.EVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.ValveTiming.IVO, self.data.linearInterpolate(self.ValveTiming.IVO, 1))
        ax.annotate('IVO %.3g $^\circ$CA' % self.ValveTiming.IVO,
                    xy=(self.ValveTiming.IVO, self.data.linearInterpolate(self.ValveTiming.IVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        # ax.scatter(self.ValveTiming.int, self.data.linearInterpolate(self.ValveTiming.int, 1))
        ax.annotate('$p_{im}$ %.4g bar' % (self.pim / 1.e5),
                    xy=(self.ValveTiming.int, self.data.linearInterpolate(self.ValveTiming.int, 1)), xycoords='data',
                    xytext=(-28, 40), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))
        ax.annotate('$p_{em}$ %.4g bar' % (self.pem / 1.e5),
                    xy=(self.ValveTiming.exh, self.data.linearInterpolate(self.ValveTiming.exh, 1)), xycoords='data',
                    xytext=(-0, 40), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.SOC, self.data.linearInterpolate(self.SOC, 1))
        ax.annotate('SOC %.3g $^\circ$CA' % self.SOC,
                    xy=(self.SOC, self.data.linearInterpolate(self.SOC, 1)), xycoords='data',
                    xytext=(-80, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))
        #
        # ax.scatter(self.EVC, table.linearInterpolate(self.EVC, 1))
        # ax.annotate('EVC %.3g $^\circ$CA' % self.EVC,
        #             xy=(self.EVC, table.linearInterpolate(self.EVC, 1)), xycoords='data',
        #             xytext=(0, 10), textcoords='offset points',
        #             arrowprops=dict(arrowstyle="->"))
        plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])

        # ax.axhline(y=0, color='r', linestyle="-.")
        ax.axvline(x=0, color='g', linestyle=":")
        ax.axvline(x=180, color='g', linestyle=":")
        ax.axvline(x=-180, color='g', linestyle=":")
        ax.axvline(x=360, color='g', linestyle=":")
        ax.axvline(x=-360, color='g', linestyle=":")

        plt.tight_layout()
        plt.show()

    def analyze(self):
        from ArrayTable import ArrayTable
        PV = ArrayTable(3, 0)
        self.data.doQuickSort(0)
        for i in range(self.data.row):
            PV.append(
                [self.data.table[0].data[i], self.CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])

        work = PV.integrate(2, _colx=1)
        print("Injected fuel {} mg".format(self.mix.gf*1.e6))

        print("BMEP={} bar".format(work / self.CylGeo.displacedVolume() / 1.e5))
        # PV.plot(1)
        # PV.plot(2, 0)
        # PV.plot(2,1)

        print("thermal efficiency {}".format(work / (self.mix.gf * 42700e3)))
