def DicplacementVolumn(Bore, Stroke, NumberOfCYlinders=1):
    import math
    return math.pi / 4. * pow(Bore, 2) * Stroke * NumberOfCYlinders

def WibeFunction(theta,SOC,TCD,a=6.908,m=2):
    from math import exp
    if theta<SOC:return 0
    else:
        return 1-exp(-a*pow((theta-SOC)/TCD,m+1))


class CylinderGeometry:
    def __init__(self, bore, stroke, connecting_rod_length=None, compression_ratio=None, number_of_cylinders=None):
        self.bore = bore
        self.stroke = stroke
        self.compression_ratio = compression_ratio
        if connecting_rod_length is None:
            print("连杆长度选择提示：")
            print("低速二冲程十字头式:{}~{}".format(3.5 * stroke / 2, 4. * stroke / 2))
            print("中速柴油机:{}~{}".format(3.8 * stroke / 2, 4.6 * stroke / 2))
            print("高速柴油机:{}~{}".format(3.5 * stroke / 2, 4.3 * stroke / 2))
        self.connecting_rod_length = connecting_rod_length
        self.__Lam = self.stroke / 2 / connecting_rod_length  # 曲柄连杆比
        if self.__Lam > 1. / 3.:
            print("请检查连杆长度是否有误，通常该连杆长度L>{}".format(3. * stroke / 2))
        print("曲柄连杆比:{}".format(self.__Lam))

        import math
        self.strokeTRC2 = math.pi * pow(self.bore, 2) / 4. * (
                self.stroke / (self.compression_ratio - 1.) + self.stroke / 2. + self.stroke / self.__Lam / 2.)
        self.strokeTRC3 = math.pi * pow(self.bore, 2) * self.stroke / 4. / 2.
        self.strokeTRC4 = self.strokeTRC3 / self.__Lam
        self.strokeTRC5 = self.strokeTRC3 * math.pi / 180
        if self.TDCclearanceHeight() < 0:
            raise Exception("TDC clearance height is less than 0!!!Parameters have some problem.")
        print("TDC clearance height is:{} mm.".format(1.e3 * self.TDCclearanceHeight()))
        print("Single cylinder displacement volume is:{}L.".format(1.e3 * self.displacedVolume()))
        if number_of_cylinders is not None:
            print("Total displacement volume of the engine is:{}L".format(
                1.e3 * number_of_cylinders * self.displacedVolume()))

    def getFi(self, V):
        """
        由容积计算曲轴转角
        :param V:单缸容积(m^3)
        :return:曲轴转角(0~180°)
        """
        if V < self.V(0) or V > self.V(180):
            raise Exception("Volume is too large to find corresponding crank angle!!")
        else:
            fi = 90
            while abs(self.V(fi) - V) / V > 1.e-5:
                fi = fi - (self.V(fi) - V) / ((self.V(fi + 0.1) - self.V(fi - 0.1)) / 2. / 0.1)
            return fi

    def V(self, Fi):
        import math
        Temp = Fi * math.pi / 180.
        H1 = math.cos(Temp)
        H2 = math.sin(Temp)
        H3 = math.sqrt(1 - pow(self.__Lam, 2) * pow(H2, 2))
        return self.strokeTRC2 - self.strokeTRC3 * H1 - self.strokeTRC4 * H3

    def DV(self, Fi):
        import math
        Temp = Fi * math.pi / 180.
        H1 = math.cos(Temp);
        H2 = math.sin(Temp);
        H3 = math.sqrt(1 - pow(self.__Lam, 2) * pow(H2, 2))
        return self.strokeTRC5 * H2 * (1 + self.__Lam * H1 / H3)

    def clearanceVolume(self):
        return self.V(0)

    def TDCclearanceHeight(self):
        import math
        return self.clearanceVolume() / (math.pi / 4 * self.bore ** 2)

    def displacedVolume(self):
        import math
        return math.pi / 4 * self.bore ** 2 * self.stroke
        # return self.V(180) - self.V(0)

    def heatExchangeArea(self, crank_angle):
        import math
        area = math.pi / 4 * self.bore ** 2
        return 2 * area + 4 * self.V(crank_angle) / self.bore

    def plotVolume(self):
        import numpy as np
        _x = np.linspace(0, 720, 1500)
        _vol = list()
        for i in _x:
            _vol.append(self.V(i) * 1.e3)
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(13, 5))
        plt.plot(_x, _vol)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Single cylinder volume(L)")
        plt.xlim((0, 720))
        # plt.ylim((self.V(0)*1.e3,self.V(180)*1.e3))
        # plt.legend()
        plt.xticks(np.linspace(0, 720, 9))
        plt.grid()
        plt.show()

    def plotHeatExchangeArea(self):
        import numpy as np
        _x = np.linspace(0, 720, 1500)
        _area = list()
        for i in _x:
            _area.append(self.heatExchangeArea(i))
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(13, 5))
        plt.plot(_x, _area)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Heat exchange area(m^2)")
        plt.xlim((0, 720))
        plt.xticks(np.linspace(0, 720, 9))
        plt.grid()
        plt.show()


class SingleWiebe:
    def __init__(self, start_of_combstion, total_combustion_duration, m=2, a=6.908):
        self.SOC = start_of_combstion
        self.TCD = total_combustion_duration
        if self.TCD <= 0:
            raise Exception("Total combustion duration can not less than 0")
        self.a = a
        self.m = m

    def X(self, crank_angle):
        if crank_angle < self.SOC: return 0
        import math
        temp = -self.a * pow(((crank_angle - self.SOC) / self.TCD), (self.m + 1.))
        return 1 - math.exp(temp)

    def DX(self, crank_angle):
        if crank_angle < self.SOC: return 0
        import math
        temp = -self.a * pow(((crank_angle - self.SOC) / self.TCD), (self.m + 1.))
        return self.a * (self.m + 1) / self.TCD * pow(((crank_angle - self.SOC) / self.TCD), (self.m + 1.)) * math.exp(
            temp)

    def plot(self):
        import numpy as np
        _crankangle = np.linspace(self.SOC - 10, self.SOC + self.TCD + 10, 200)
        y1 = list();
        y2 = list()
        for i in _crankangle:
            y1.append(self.DX(i))
            y2.append(self.X(i))

        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(_crankangle, y1)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()

        plt.subplot(122)
        plt.plot(_crankangle, y2)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()


class DoubleWiebe:
    def __init__(self, PF, PSOC, PTCD, DTCD, PM=2, DM=0.8, DSOCin=None):
        if PF > 1 or PF < 0:
            raise Exception("Premixed fraction must be less than 1 and greater than 0, but {} was given".format(PF))
        self.PF = PF
        self.DF = 1 - self.PF
        self.PSOC = PSOC
        self.PTCD = PTCD
        self.PM = PM
        self.DTCD = DTCD
        self.DM = DM
        if DSOCin is None:
            self.DSOC = self.PSOC + self.PTCD / 2.
        else:
            self.DSOC = DSOCin
        self.PremixWiebe = SingleWiebe(self.PSOC, self.PTCD, self.PM)
        self.DisfusionWiebe = SingleWiebe(self.DSOC, self.DTCD, self.DM)

    def X(self, CA):
        return self.PF * self.PremixWiebe.X(CA) + self.DF * self.DisfusionWiebe.X(CA)

    def DX(self, CA):
        return self.PF * self.PremixWiebe.DX(CA) + self.DF * self.DisfusionWiebe.DX(CA)

    def getData(self):
        import ArrayTable as AT
        import numpy as np
        data = AT.ArrayTable(3, 0)
        _crankangle = np.linspace(self.PSOC - 10, self.DSOC + self.DTCD + 20, 200)
        for i in _crankangle:
            data.append([i, self.DX(i), self.X(i)])
        return data

    def plot(self):
        import numpy as np
        _crankangle = np.linspace(self.PSOC - 10, self.DSOC + self.DTCD + 10, 200)
        y1 = list();
        y2 = list()
        for i in _crankangle:
            y1.append(self.DX(i))
            y2.append(self.X(i))

        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(_crankangle, y1)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()

        plt.subplot(122)
        plt.plot(_crankangle, y2)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()


class HeatReleaseData:
    def __init__(self, _data):
        self.data = _data
        self.SOC = self.__findSOC()
        print("Start of combustion:{}".format(self.SOC))
        self.EOC = self.__findEOC()
        print("End of combustion:{}".format(self.EOC))
        self.__analysis()
        self.SOCIndex = self.data.table[0].data.index(self.SOC)
        self.EOCIndex = self.data.table[0].data.index(self.EOC)

    def __findSOC(self):
        index = self.data.findMaxValueIndex(1)
        try:
            while (self.data.table[1].data[index] > 1.e-8):
                index -= 1;
        except:
            print("Can not find start of combustion, will return the first element")
            return self.data.table[0].data[-1]
        return self.data.table[0].data[index]

    def __findEOC(self):
        index = self.data.findMaxValueIndex(1)
        try:
            while (self.data.table[1].data[index] > 1.e-4):
                index += 1;
        except:
            print("Can not find end of combustion, will return the last element")
            return self.data.table[0].data[-1]
        return self.data.table[0].data[index]

    def __analysis(self):
        from ArrayTable import PhsicalVarList
        burnedfraction = PhsicalVarList([], "Total burned fraction", "/")
        self.CA10 = self.CA50 = self.CA90 = None
        print("Maximum burned fraction is {}".format(self.data.integrate(1)))
        index = 0

        while True:
            burnedfraction.data.append(self.data.integrate(1, index, 0))
            index += 1
            if index >= self.data.row:
                print("Heat release data has some problem, can not find CA10,CA50,CA90")
                return

        self.CA10 = self.data.table[0].data[index]
        print("Find CA10 at {} °CA".format(self.CA10))

        while self.data.integrate(1, index, 0) < 0.5:
            index += 1
            if index >= self.data.row:
                self.data.appendColumn(burnedfraction)
                print("Heat release data has some problem, can not find CA50,CA90")
                return
        self.CA50 = self.data.table[0].data[index]
        print("Find CA50 at {} °CA".format(self.CA50))

        while self.data.integrate(1, index, 0) < 0.9:
            index += 1
            if index >= self.data.row:
                print("Heat release data has some problem, can not find CA90")
                return
        self.CA90 = self.data.table[0].data[index]
        print("Find CA90 at {} °CA".format(self.CA90))

    def plot(self):
        import matplotlib.pyplot as plt
        fig1, ax1 = plt.subplots(211, figsize=(10, 10))

        fig1, ax1 = plt.subplots(211, figsize=(10, 10))

    def regressWithSingleWiebe(self):
        from sko.PSO import PSO
        pso = PSO(func=self.__regressWithSingleWiebeFunc, dim=1, pop=100, max_iter=150, lb=[0], ub=[20])
        fitness = pso.run()
        print("exp(-{}(x-{})/{})^{}".format(6.908, self.SOC, self.EOC - self.SOC, pso.gbest_x[0]))
        print('Absolute error=', pso.gbest_y[0])
        from ArrayTable import ArrayTable
        from Cylinder import SingleWiebe
        result = ArrayTable(3, 0)
        wiebe = SingleWiebe(self.SOC, self.EOC - self.SOC, pso.gbest_x[0])
        for i in range(0, self.data.row):
            result.append(
                [self.data.table[0].data[i], self.data.table[1].data[i], wiebe.DX(self.data.table[0].data[i])])
        result.plot(0, [1, 2])
        return result

    def regressWithDoubleWiebe(self):
        from sko.PSO import PSO
        pso = PSO(func=self.__regressWithDoubleWiebeFunc, dim=5, pop=100, max_iter=400, lb=[0.0, 10, 10, 0.1, 0.1],
                  ub=[1, 180, 180, 5, 5], w=0.8, c1=2, c2=2)
        fitness = pso.run()
        print("{}*exp(-6.908(x-{})/{})^{}+{}*exp(-6.908(x-{})/{})^{})".format(pso.gbest_x[0], self.SOC, pso.gbest_x[1],
                                                                              pso.gbest_x[3], 1 - pso.gbest_x[0],
                                                                              self.EOC - pso.gbest_x[2], pso.gbest_x[2],
                                                                              pso.gbest_x[4]))
        print('best_x is ', pso.gbest_x)
        print('best_y is ', pso.gbest_y)
        from ArrayTable import ArrayTable
        from Cylinder import DoubleWiebe
        result = ArrayTable(3, 0)
        wiebe = DoubleWiebe(pso.gbest_x[0], self.SOC, pso.gbest_x[1], pso.gbest_x[2], pso.gbest_x[3], pso.gbest_x[4],
                            self.EOC - pso.gbest_x[2])
        for i in range(0, self.data.row):
            result.append(
                [self.data.table[0].data[i], self.data.table[1].data[i], wiebe.DX(self.data.table[0].data[i])])
        result.setTableHeader(["Crank angle", "Original HRR", "PSO HRR"])
        result.setTableUnit(["CA", "1/deg", "1/deg"])
        result.plot(0, [1, 2])
        return result

    def __regressWithSingleWiebeFunc(self, m, a=6.908):
        from Cylinder import SingleWiebe
        wiebe = SingleWiebe(self.SOC, self.EOC - self.SOC, m, a)
        result = 0
        for i in range(self.SOCIndex, self.EOCIndex):
            result += abs(wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i])
        return result

    def __regressWithDoubleWiebeFunc(self, PF, PTCD, DTCD, PM=2, DM=0.8):
        from Cylinder import DoubleWiebe
        wiebe = DoubleWiebe(PF, self.SOC, PTCD, DTCD, PM, DM, self.EOC - DTCD)
        result = 0
        for i in range(self.SOCIndex, self.EOCIndex):
            result += abs(wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i])
        return result


class CylinderPressure:
    def __init__(self, _data):
        self.data = _data

    def FFTTrasform(self, n, tau=4.):
        """
        傅里叶变换分析
        :param n:
        :param tau:
        :return:
        """

        def basicFrequency(n, CAInterval, stroke=4.):
            t_cycle = 30. * stroke / float(n)
            t_interval = CAInterval * t_cycle / 720.
            return 1 / t_interval

        CAInterval = self.data.table[0].data[1] - self.data.table[0].data[0]
        baseF = basicFrequency(n, CAInterval, tau)
        print("Data collection interval is {} CA.".format(CAInterval))
        print("Basic frequency is {} kHz".format(baseF / 1.e3))
        from numpy import fft
        import numpy as np
        dofft = fft.fft(self.data.table[1].data)
        from ArrayTable import ArrayTable
        result = ArrayTable(3, 0)
        result.setTableHeader(["Frequency", "Amplitude", "Phase"])
        result.setTableUnit(["Hz", "/", ""])
        for i in range(self.data.row):
            result.append([i * baseF, np.abs(dofft[i]), np.angle(dofft[i])])
        return result

    def FFTFilter(self, _cutOffFrequency, amplitudeMutiplier=1):
        from numpy import fft
        from ArrayTable import ArrayTable
        import numpy as np
        dofft = fft.fft(self.data.table[1].data)
        temp = dofft
        for i in range(len(dofft)):
            if i > _cutOffFrequency:
                temp[i] = 0
        temp = [each * amplitudeMutiplier for each in temp]
        result = ArrayTable(3, 0)
        result.table[1].data = np.abs(fft.ifft(temp))
        result.table[0].data = self.data.table[0].data
        result.table[2].data = self.data.table[1].data
        return result

    def PVDiagram(self, CylGeo):
        from ArrayTable import ArrayTable
        result = ArrayTable(2, 0)
        for i in range(self.data.row):
            result.append([CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])
        result.setTableHeader(["Cylinder volume", "Cylinder pressure"])
        result.setTableUnit(["$m^3$", "bar"])
        return result

    def Temperature(self, CylGeo):
        from ArrayTable import ArrayTable
        from GasProperty import Rg
        result = ArrayTable(4, 0)
        for i in range(self.data.row):
            V = CylGeo.V(self.data.table[0].data[i])
            T = self.data.table[1].data[i] * 1.e5 * V / Rg()
            result.append(
                [self.data.table[0].data[i], T, CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])
        result.setTableHeader(["Crank angle", "Cylinder temperature", "Cylinder volume", "Cylinder pressure"])
        result.setTableUnit(["$^\circ CA$", "K", "$m^3$", "bar"])
        return result

    def startOfCombustion(self, type=0):
        self.slice(-60, 60)
        self.data.diff(1)
        self.data.diff(2)
        if type == 0:
            zeroindex = self.data.findMaxValueIndex(1)
            while self.data.table[0].data[zeroindex] > 0:
                zeroindex -= 1
            while self.data.table[3].data[zeroindex] > 0:
                zeroindex -= 1
            fi0 = self.data.table[0].data[zeroindex]
            fi1 = self.data.table[0].data[zeroindex + 1]
            p20 = self.data.table[3].data[zeroindex]
            p21 = self.data.table[3].data[zeroindex + 1]
            soc = fi0 - p20 * (fi1 - fi0) / (p21 - p20)
            print("Find start of combustion at {} CA".format(soc))
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.plot(self.data.table[0].data, self.data.table[3].data)
            plt.xlabel(self.data.table[0].ColName + "(" + self.data.table[0].ColUnit + ")")
            plt.ylabel(self.data.table[3].ColName + "(" + self.data.table[3].ColUnit + ")")
            ax.axhline(y=0, color='r', linestyle="-.")
            ax.axvline(x=0, color='g', linestyle=":")
            ax.scatter(soc, 0, color="r")
            ax.annotate('Start of combustion %.3g $^\circ$CA' % soc,
                        xy=(soc, 0), xycoords='data',
                        xytext=(-170, 50), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))
            plt.tight_layout()
            plt.show()
            return soc

        if type == 1:
            index = self.data.findMaximumDataIndex(3, 2)
            trueindex = [i for i in index if self.data.table[0].data[i] < 0]

            print("Premixed combustion starts at {} CA".format(self.data.table[0].data[trueindex[-2]]))
            print("Diffusion combustion starts at {} CA".format(self.data.table[0].data[trueindex[-1]]))

            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.plot(self.data.table[0].data, self.data.table[3].data)
            plt.xlabel(self.data.table[0].ColName + "(" + self.data.table[0].ColUnit + ")")
            plt.ylabel(self.data.table[3].ColName + "(" + self.data.table[3].ColUnit + ")")
            ax.axhline(y=0, color='r', linestyle="-.")
            ax.axvline(x=0, color='g', linestyle=":")
            ax.scatter(self.data.table[0].data[trueindex[-1]], self.data.table[3].data[trueindex[-1]], color="r")
            ax.scatter(self.data.table[0].data[trueindex[-2]], self.data.table[3].data[trueindex[-2]], color="r")

            ax.annotate('Start of combustion %.3g $^\circ$CA' % self.data.table[0].data[trueindex[-2]],
                        xy=(self.data.table[0].data[trueindex[-2]], self.data.table[3].data[trueindex[-2]]),
                        xycoords='data',
                        xytext=(-200, 0), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.annotate('Start of combustion %.3g $^\circ$CA' % self.data.table[0].data[trueindex[-1]],
                        xy=(self.data.table[0].data[trueindex[-1]], self.data.table[3].data[trueindex[-1]]),
                        xycoords='data',
                        xytext=(20, 50), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            plt.tight_layout()
            plt.show()
            return trueindex[-2], trueindex[-1]

        if type == 2:
            self.data.diff(3)
            index = self.data.findMaximumDataIndex(4, 2)
            trueindex = [i for i in index if self.data.table[0].data[i] < 0]

            print("Premixed combustion starts at {} CA".format(self.data.table[0].data[trueindex[-2]]))
            print("Diffusion combustion starts at {} CA".format(self.data.table[0].data[trueindex[-1]]))

            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.plot(self.data.table[0].data, self.data.table[4].data)
            plt.xlabel(self.data.table[0].ColName + "(" + self.data.table[0].ColUnit + ")")
            plt.ylabel(self.data.table[4].ColName + "(" + self.data.table[4].ColUnit + ")")
            ax.axhline(y=0, color='r', linestyle="-.")
            ax.axvline(x=0, color='g', linestyle=":")
            ax.scatter(self.data.table[0].data[trueindex[-1]], self.data.table[4].data[trueindex[-1]], color="r")
            ax.scatter(self.data.table[0].data[trueindex[-2]], self.data.table[4].data[trueindex[-2]], color="r")

            ax.annotate('Start of combustion %.3g $^\circ$CA' % self.data.table[0].data[trueindex[-2]],
                        xy=(self.data.table[0].data[trueindex[-2]], self.data.table[4].data[trueindex[-2]]),
                        xycoords='data',
                        xytext=(-200, 0), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))
            ax.annotate('Start of combustion %.3g $^\circ$CA' % self.data.table[0].data[trueindex[-1]],
                        xy=(self.data.table[0].data[trueindex[-1]], self.data.table[4].data[trueindex[-1]]),
                        xycoords='data',
                        xytext=(20, 50), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))
            plt.tight_layout()
            plt.show()
            return trueindex[-2], trueindex[-1]

    def slice(self, left=None, right=None):
        if left is None or left < self.data.table[0].data[0]:
            left = self.data.table[0].data[0]
        if right is None or right > self.data.table[0].data[-1]:
            right = self.data.table[0].data[-1]

        # maxdataindex=self.data.findMaxValueIndex(1)
        while self.data.table[0].data[-1] > right:
            self.data.popRow()

        while self.data.table[0].data[0] < left:
            self.data.delRow(0)

    def sliceToComAndExhaust(self, IVC=None, EVO=None):
        maxPreIndex = self.data.findMaxValueIndex(1)
        if IVC is None:
            IVC = self.data.table[0].data[maxPreIndex] - 180. + 5.
        if EVO is None:
            EVO = self.data.table[0].data[maxPreIndex] + 180. - 5.
        self.slice(IVC, EVO)

    def ployTropicIndex(self, CylGeo):
        from ArrayTable import ArrayTable
        from math import log
        result = ArrayTable(4, 0)
        result.setTableHeader(["Crank angle", "poly tropic index", "V", "p"])
        result.setTableUnit(["°CA", "/", "", ""])
        for i in range(self.data.row - 1):
            temp = log(self.data.table[1].data[i + 1] / self.data.table[1].data[i]) / log(
                CylGeo.V(self.data.table[0].data[i]) / CylGeo.V(self.data.table[0].data[i + 1]))
            result.append([self.data.table[0].data[i], temp,
                           CylGeo.V(self.data.table[0].data[i]) / CylGeo.V(self.data.table[0].data[i + 1]),
                           self.data.table[1].data[i + 1] / self.data.table[1].data[i]])

        return result

    def netHeatReleaseRate(self):
        def gamma(T):
            return 1.338 - 6.0e-5 * T + 1.0e-8 * T * T

    def smooth(self, smoothType=None):
        smoothTypelist = [None, "five points three times smooth", "FFT smooth"]
        if smoothType not in smoothTypelist:
            print("Invalid smooth type!!")
            print("allowed smooth type are: ")
            print(smoothTypelist)
            raise Exception("")

        def fivePointsThreeTimesfun(data):
            result = [0] * len(data)
            i = 0;
            result[i] = (69. * data[i] + 4. * data[i + 1] - 6. * data[i + 2] + 4. * data[i + 3] - data[i + 4]) / 70.
            i = 1;
            result[i] = (2. * data[i - 1] + 27. * data[i] + 12. * data[i + 1] - 8. * data[i + 2] + 2. * data[
                i + 3]) / 35.
            i = -2;
            result[i] = (2. * (data[i - 3] + data[i + 1]) - 8. * data[i - 2] + 12. * data[i - 1] + 27. * data[i]) / 35.
            i = -1;
            result[i] = (-data[i - 4] + 4. * (data[i - 3] + data[i - 1]) - 6. * data[i - 2] + 69. * data[i]) / 70.
            for i in range(2, len(data) - 2):
                result[i] = (-3. * (data[i - 2] + data[i + 2]) + 12. * (data[i - 1] + data[i + 1]) + 17. * data[
                    i]) / 35.
            return result

        if smoothType == "five points three times smooth":
            aftersmooth = fivePointsThreeTimesfun(self.data.table[1].data)
            from ArrayTable import PhsicalVarList
            aftersmoothcol = PhsicalVarList()
            aftersmoothcol.ColName = "Pressure after smooth"
            aftersmoothcol.ColUnit = self.data.table[1].ColUnit
            aftersmoothcol.data = aftersmooth
            self.data.table.append(aftersmoothcol)
            self.data.col += 1
            return self.data


from Compressor import *
from GasProperty import *


class MillerCycle:
    def __init__(self, CylGeo, p0=101325, T0=273.15 + 20):
        self.p0 = p0
        self.T0 = T0
        self.CylGeo = CylGeo
        from ArrayTable import ArrayTable
        table = ArrayTable(5, 0)
        table.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        table.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])
        self.data = table

    def initCompressor(self, etak, pik):
        self.pik = pik
        self.etak = etak
        self.pk = self.p0 * pik
        self.Tk = TAfterCompressor(self.T0, pik, etak)

    def intake(self, IVC=-30):
        # 这里IVC为进气门相对于下止点的关闭角度，为0时代表在下止点关闭
        self.IVC = IVC + 180
        print("Effective compression ratio is {}".format(self.CylGeo.V(180 + IVC) / self.CylGeo.V(0)))
        if IVC <= -90 or IVC > 150:
            raise Exception("Intake valve close timing can not earlier than -90 and later than 150 with respect to "
                            "BDC,but {} is given".format(IVC))
        from GasProperty import Rg
        from numpy import arange
        for i in arange(0, 180 + IVC, 1):
            V = self.CylGeo.V(i)
            self.data.append([i, self.CylGeo.V(i), self.pk, self.Tk, self.pk * V / Rg(1.e8) / self.Tk])

    def adiabaticompress(self):
        from numpy import arange
        for i in arange(self.IVC + 0.1, 360., 1):
            k = k_Justi(self.data.table[3].data[-1])
            V = self.CylGeo.V(i)
            p = pow(self.data.table[1].data[-1] / V, k) * self.data.table[2].data[-1]
            T = pow(self.data.table[1].data[-1] / V, k - 1) * self.data.table[3].data[-1]
            self.data.append([i, V, p, T, self.data.table[4].data[-1]])

    def preBurn(self, Rp, alpha=1, Hu=42700e3, L0=14.3):
        if Rp < 0 or Rp > 1: raise Exception("premixed burn fraction value error!!")
        self.Rp = Rp
        self.L0 = L0
        self.alpha = alpha
        from GasProperty import DieselMixture
        mix = DieselMixture()
        mix.init_With_Ma_r(self.data.table[4].data[-1], 0, self.data.table[4].data[-1] / alpha / L0)
        self.mix = mix

        def fun(x, T):
            return (Hu - mix.u(x, T)) / (alpha * L0 + x) / mix.cv(x, T)

        x = 0
        self.m3 = self.data.table[4].data[-1]
        T = self.data.table[3].data[-1]
        step = Rp / 50.
        while x < Rp:
            T += fun(x, T) * step
            x += step
            m = self.m3 + x * self.m3 / alpha / L0
            p = m * mix.Rg_gas(x) * T / self.CylGeo.V(360.)
            self.data.append([360., self.CylGeo.V(360.), p, T, m])

    def premixBurn(self, Rp, alpha=1, Hu=42700e3, L0=14.3):
        self.Rp = Rp
        self.L0 = L0
        self.alpha = alpha
        from GasProperty import DieselMixture
        mix = DieselMixture()
        mix.init_With_Ma_r(self.data.table[4].data[-1], 0, self.data.table[4].data[-1] / alpha / L0)
        self.mix = mix

        T3 = self.data.table[3].data[-1]
        p3 = self.data.table[2].data[-1]
        T4 = T3 + 100

        def fun(Tz):
            return (cv_mean_Justi(T3, Tz, mix.AFAK(Rp)) * (Rp + L0 * alpha) * (Tz - T3) - Rp * Hu) / 1.e3

        h = 1.
        while abs(fun(T4)) > 1:
            dfun = (fun(T4 + h) - fun(T4 - h)) / 2. / h
            T4 -= fun(T4) / dfun
        self.T4 = T4
        print("T4={}".format(T4))
        self.p4 = p3 * T4 / T3
        self.m4 = self.p4 * self.CylGeo.V(360.) / Rg(mix.AFAK(Rp)) / T4
        self.data.append([360., self.CylGeo.V(360.), self.p4, self.T4, self.m4])

    def DisffusionBurn(self, Hu=42700e3):
        def fun(T5):
            return ((1 + self.L0 * self.alpha) * cp_mean_Justi(self.T4, T5) * (T5 - self.T4) + (
                    self.Rp + self.alpha * self.L0) * Rg(self.alpha) * (T5 - self.T4) - (1 - self.Rp) * Hu) / 1.e4

        T5 = self.T4 + 100
        h = 1.
        while abs(fun(T5)) > 1:
            dfun = (fun(T5 + h) - fun(T5 - h)) / 2. / h
            T5 -= fun(T5) / dfun
        self.T5 = T5
        print("T5={}".format(T5))
        self.V5 = self.T5 / self.T4 * self.CylGeo.V(0)
        self.data.append([self.V5, self.p4, self.T5, self.p4 * self.V5 / Rg(self.alpha) / self.T5])

    def disffBurn(self, Hu=42700e3):
        def fun(x, T):
            return (Hu - self.mix.h(x, T)) / (self.mix.cp(x, T) * (self.alpha * self.L0 + x))

        x = self.Rp
        T = self.data.table[3].data[-1]
        T4 = self.data.table[3].data[-1]
        V4 = self.CylGeo.V(360.)
        p4 = self.data.table[2].data[-1]
        step = (1 - self.Rp) / 100.
        while x < 1:
            T += fun(x, T) * step
            x += step
            m = self.m3 + x * self.m3 / self.alpha / self.L0
            V = T / T4 * V4
            self.data.append([self.CylGeo.getFi(V) + 360., V, p4, T, m])

    def expansion(self):
        from numpy import arange
        for i in arange(self.data.table[0].data[-1], 540., 1):
            k = k_Justi(self.data.table[3].data[-1])
            V = self.CylGeo.V(i)
            p = pow(self.data.table[1].data[-1] / V, k) * self.data.table[2].data[-1]
            T = pow(self.data.table[1].data[-1] / V, k - 1) * self.data.table[3].data[-1]
            self.data.append([i, V, p, T, self.data.table[4].data[-1]])

        self.p6 = self.data.table[2].data[-1]
        self.T6 = self.data.table[3].data[-1]

    def supercriticalExhaust(self, etat):
        from Compressor import piT
        piTinit = 2.

        def fun(piTx):
            return piTx - piT(self.pik, self.etak * etat, self.p0 * piTx / self.p6 * self.T6, self.T0,
                              self.alpha * self.L0)
            # return piK(self.p0 * piT / self.p6 * self.T6, piT, etat, self.etak, self.T0,
            #            self.alpha * self.L0) - self.pik

        h = 0.01
        while abs(fun(piTinit)) > 1.e-5:
            dfun = (fun(piTinit + h) - fun(piTinit - h)) / 2 / h
            piTinit -= fun(piTinit) / dfun
            print(piTinit)

        # h = 0.01
        # while abs(fun(piTinit)) > 1.e-2:
        #     dfun = (fun(piTinit + h) - fun(piTinit - h)) / 2 / h
        #     piTinit -= fun(piTinit) / dfun
        #     print(piTinit)
        self.T7 = self.p0 * piTinit / self.p6 * self.T6
        self.p7 = self.p0 * piTinit
        self.data.append(
            [540., self.CylGeo.V(540), self.p7, self.T7, self.p7 * self.CylGeo.V(540) / Rg(self.alpha) / self.T7])

    def subcritical(self):
        from numpy import arange
        for i in arange(540 + 0.1, 720, 1):
            V = self.CylGeo.V(i)
            self.data.append([i, V, self.p7, self.T7, self.p7 * V / Rg(self.alpha) / self.T7])


def IdealMillerCylcle(T2=400, Rp=0.01, Hu=42700e3, alpha=1.1, L0=14.3, k=None, eps=18, epsm=15, pik=4, etaTK=0.5,p0=1.e5):
    from GasProperty import cv_Justi, cp_Justi, k_Justi
    if k is None:
        k=k_Justi(T2)
    cv = cv_Justi(T2)
    lamda = Rp * Hu / (alpha * L0 * cv * pow(epsm, k - 1) * T2) + 1
    cp = cp_Justi(T2)
    rho = Hu * (1 - Rp) / (k * Hu * Rp + alpha * L0 * pow(epsm, k - 1) * cp * T2) + 1
    pit = pow(1 - pow(eps / epsm, k - 1) * (pow(pik, (k - 1) / k) - 1) / lamda / pow(rho, k) / etaTK, -k / (k - 1))
    eta = 1 - (lamda * pow(rho, k) / pow(eps, k - 1) + (k - 1) * eps / pow(epsm, k) - k / pow(epsm, k - 1) - (k - 1) * (
                1 - pit / pik) * (eps - 1) / pow(epsm, k)) / (lamda - 1 + k * lamda * (rho - 1))
    # BMEP=p0*pik*(((lamda-1)+k*lamda*(rho-1))
    return eta


def x(theta, start=0, end=60):
    from math import cos, pi
    return 0.5 * (1 - cos(pi * (theta - start) / (end - start)))


def FMEP(D, cm, pme) -> float:
    """
    计算平均机械损失压力，适用于四冲程
    :param D:Bore，气缸直径(m)
    :param cm:Mean piston moving velocity，活塞平均运动速度(m/s)
    :param pme:Brake mean effective pressure，平均有效压力(Pa)
    :rtype: FMEP，平均机械损失压力(Pa)
    """
    return pow(D, -0.1778) * (0.00855 * cm + 0.0789 * pme / 1.e6 - 0.0214) * 1.e6


def MeEff(D, cm, pme) -> float:
    """
    计算机械效率，适用于四冲程
    :param D:Bore，气缸直径(m)
    :param cm:Mean piston moving velocity，活塞平均运动速度(m/s)
    :param pme:Brake mean effective pressure，平均有效压力(Pa)
    :rtype: Mechanical efficiency,机械效率=BMEP/(BMEP+FMEP)(Pa)
    """
    fmep = FMEP(D, cm, pme)
    return pme / (pme + fmep)


def BMEP(n, Pe, TVH, stroke=4):
    """
    计算平均有效压力
    :param n: Speed,发动机转速(r/min)
    :param Pe: Brake power 发动机功率(W)
    :param TVH: Total displacement volume,发动机总排量(m^3)
    :param stroke: 发动机冲程
    :return: Brake mean effective pressure,平均有效压力(Pa)
    """
    return Pe * (30 * stroke / n) / TVH


def BSFC(etai, etam, Hu=42496e3):
    """
    计算燃油消耗率
    :param etai:指示热效率
    :param etam: 机械效率
    :param Hu: 燃料低热值(J/kg)
    :return: 有效燃油消耗率
    """
    print("低速柴油机165~190g/(kW·h)")
    print("中速柴油机195~220g/(kW·h)")
    print("高速柴油机210~230g/(kW·h)")
    return 3.6e6 / (Hu / 1.e3 * etai * etam)


def BSFCexample(n):
    def fun(x):
        return -827.92 * x ** 4 + 2163.5 * x ** 3 - 1835.2 * x ** 2 + 526.25 * x + 177.48

    from ArrayTable import ArrayTable
    result = ArrayTable(2, 0)
    result.setTableHeader(["Speed", "BSFC"])
    result.setTableUnit(["rpm", "g/(kW*h)"])
    from numpy import arange
    for i in arange(0.3, 1, 0.01):
        result.append([n * i, fun(i)])
    return result


def mfcycle(ge, Pe, n, i=1):
    """
    计算单缸循环喷油量
    :param ge: BSFC，燃油消耗率(g/(kW*h))
    :param Pe: Brake power 发动机功率(kW)
    :param n: Speed，发动机转速(r/min)
    :param i: Number of cylinders，气缸数，i=1时计算循环的喷油量
    :return: fuel mass per cylinder per cycle,单缸循环喷油量(kg)
    """
    result = ge * Pe / 3600 * 120 / n / i * 1.e-3
    print("fuel mass injected per cylinder per cycle:{}mg".format(result * 1.e6))
    return result

def mfcylclePim(VE):
    pass

def massAir(Vd,ps=1.e5,Ts=300,VE=1.0):
    from GasProperty import Rg
    return VE*ps*Vd/Rg()/Ts

def thermalUsageIndex(gi, alpha, Cm, D, S, n, Ps, Ts):
    """
    热利用系数
    :param gi:ISFC,指示燃油消耗率(g/(kW*h))
    :param alpha:Excess air fuel ratio,过量空气系数
    :param Cm:mean piston velocity,活塞平均运动速度(m/s)
    :param D:Bore,缸径(m)
    :param S:Stroke,冲程(m)
    :param n:Speed,发动机转速(r/min)
    :param Ps:Intake manifold pressure(Pa),进气管压力
    :param Ts:Intake manifold temperature(K),进气管温度
    :return:热利用系数
    """
    R = pow(alpha * gi / 1e3, 0.5) * pow(Ts, 1.5) * pow(Cm, 0.78) * (0.5 + D / (2 * S)) / (
            (D * n) * pow(D * Ps / 1.e6, 0.22))
    return 1.028 - 0.00096 * R


def simpleWoschini():
    pass
