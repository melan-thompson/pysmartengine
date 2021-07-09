from ArrayTable import ArrayTable
import math

def flowUnitArea(p1, T1, R1, k1, p2, T2=0, R2=0, k2=0):
    if p1 < 0 or p2 < 0:
        raise Exception("Pressure can not less than 0!!!")
    if T1 < 0 or T2 < 0:
        raise Exception("Temperature can not less than 0!!!")
    if R1 < 0 or R2 < 0:
        raise Exception("Gas constant can not less than 0!!!")
    if k1 < 0 or k2 < 0:
        raise Exception("Inappropriate specific heat ratio!!!")
    from math import sqrt, pow
    TMP2 = p2 / p1
    if abs(TMP2 - 1) < 1.e-2:
        return 0
    elif TMP2 < 1.0:  ##正向流动
        TMP1 = sqrt(2 * k1 / R1 / T1)
        TMP3 = pow(2 / (k1 + 1), 1 / (k1 - 1))
        if TMP2 > pow(TMP3, k1):  ##正向亚临界流动
            return TMP1 * p1 * sqrt((pow(TMP2, 2 / k1) - pow(TMP2, (k1 + 1) / k1)) / (k1 - 1))
        else:  # 正向超临界流动
            return TMP1 * p1 / sqrt(k1 + 1) * TMP3
    else:  # 反向流动
        TMP1 = sqrt(2 * k2 / R2 / T2)
        TMP2 = p1 / p2
        TMP3 = pow(2 / (k2 + 1), 1 / (k2 - 1))
        if TMP2 > pow(TMP3, k2):  # 反向亚临界流动
            return -TMP1 * p2 * sqrt((pow(TMP2, 2 / k2) - pow(TMP2, (k2 + 1) / k2)) / (k2 - 1))
        else:  # 反向超临界流动
            return -TMP1 * p2 / sqrt(k2 + 1) * TMP3


def mod(angle, TDC=0):
    left = TDC - 360
    right = TDC + 360
    if left <= angle < right:
        return angle
    while angle < left:
        angle += 720
    while angle >= right:
        angle -= 720
    return angle


class ValveSimple:
    def __init__(self, _flowArea=1.):
        from ArrayTable import ArrayTable
        self.flowArea = _flowArea
        self.data = ArrayTable(3, 0)
        self.data.setTableHeader(["time", "Pressure ratio", "Mass flow rate"])
        self.data.setTableUnit(["s", "/", "kg/s"])

        # 连接关系
        self.next = None
        self.last = None

    def __massFlowRate(self, p1, T1, R1, k1, p2, T2=0, R2=0, k2=0):
        return self.flowArea * flowUnitArea(p1, T1, R1, k1, p2, T2, R2, k2)

    def connect_to(self, last, next):
        self.last = last
        self.last.next = self

        self.next = next
        self.next.last = self

    def update(self, t):
        self.dm_record = self.__massFlowRate(self.last.p, self.last.T, self.last.Rg, self.last.k, self.next.p,
                                             self.next.T, self.next.Rg, self.next.k)
        self.record(t)
        return self.dm_record

    def record(self, t):
        self.data.append([t, self.next.p / self.last.p, self.dm_record])

    def airFlowExample(self, p1, T1, p2=None, T2=None, K2=None):
        from GasProperty import Rg, k_Justi
        R1 = Rg(T1)
        R2 = Rg(T2)
        k1 = k_Justi(T1)
        k2 = k_Justi(T2)
        step = 1.e-3 * p1
        import numpy as np
        for i in np.arange(0, p1, step):
            self.data.append([i / p1, self.__massFlowRate(p1, T1, R1, k1, i)])

        if p2 is not None and T2 is not None:
            step2 = 0.01 * p2
            for i in np.arange(0.5 * p2, p2, step):
                self.data.append([p2 / i, self.__massFlowRate(i, T1, R1, k1, p2, T2, R2, K2)])

        return self.data


class Valve:
    def __init__(self, valveDiameter, open_timing_angle, close_timing_angle, _valve_lift, sigma=0, rock=1,
                 _valve_flow_coeff=None):
        self.valve_diameter = valveDiameter
        self.open_timing_angle = open_timing_angle
        self.close_timing_angle = close_timing_angle
        self.sigma = sigma * math.pi / 180.

        ## 放缩
        open = _valve_lift.table[0].data[0]
        close = _valve_lift.table[0].data[-1]
        temp = (close_timing_angle - open_timing_angle) / (close - open)
        self.valve_lift = ArrayTable(2, 0)
        self.valve_lift.setTableUnit(["CA", "m"])
        self.valve_lift.setTableHeader(["Crank angle", "Valve lift"])
        for i in range(_valve_lift.row):
            self.valve_lift.append(
                [mod(self.open_timing_angle + temp * (_valve_lift.table[0].data[i] - _valve_lift.table[0].data[0])),
                 _valve_lift.table[1].data[i] * rock])

        self.valve_lift.doQuickSort()

        # 流量系数
        self.flow_coeff = _valve_flow_coeff
        if _valve_flow_coeff is None:
            self.flow_coeff = ArrayTable(2, 0)
            self.flow_coeff.setTableHeader(["Valve lift", "Flow coefficient"])
            self.flow_coeff.setTableUnit(["m", "/"])
            maxlift = self.valve_lift.findMaxValue(1)
            print("maximum valve lift is {}".format(maxlift))
            if valveDiameter<maxlift/0.544:
                raise Exception("valve diameter is too small, it can not be less than {}".format(maxlift/0.544))
            import numpy as np
            for i in np.arange(0, maxlift, maxlift/100):
                flowcoeff=self.__flowcoefficient(i)
                if flowcoeff < 0:
                    raise Exception("valve flow coefficient can not less than 0")
                self.flow_coeff.append([i, flowcoeff])

        # 有效流通面积
        self.flowAreaData=ArrayTable(2,0)
        self.flowAreaData.setTableHeader(["Valve lift", "Flow area"])
        self.flowAreaData.setTableUnit(["$^\circ$CA", "$m^2$"])
        for i in np.arange(-540,540):
            self.flowAreaData.append([i,self.flowArea(i)])

        # 下游上游
        self.last = None
        self.next = None

        self.dm_record=0

    # valvelift从0开始
    def valveLiftStartFromZero(self):
        temp = self.valve_lift
        startpoint = self.valve_lift.table[0].data[0]
        temp.table[0] + (-startpoint)
        return temp

    # 推荐的流量系数
    def __flowcoefficient(self, _valvelift):
        return 0.98 - 3.3 * pow((_valvelift / self.valve_diameter), 2)

    # 计算流通面积
    def flowArea(self, Fi):
        from math import pi, cos, sin
        hv = self.valve_lift.linearInterpolate(Fi, 1,1)
        flowcoeff = self.__flowcoefficient(hv)
        if flowcoeff<0:
            raise Exception("valve flow coefficient can not less than 0")
        return flowcoeff * pi * hv * cos(self.sigma) * (self.valve_diameter + hv * sin(self.sigma) * cos(self.sigma))

    # 计算流量kg/s
    def massFlowRate(self, Fi,speed=None):
        # 流量单位为kg/s
        if speed is None:
            self.dm_record = flowUnitArea(self.last.p, self.last.T, self.last.Rg, self.last.k, self.next.p,
                                          self.next.T, self.next.Rg, self.next.k)*self.flowArea(Fi)
        else:
            self.dm_record = flowUnitArea(self.last.p, self.last.T, self.last.Rg, self.last.k, self.next.p,
                                          self.next.T, self.next.Rg, self.next.k) * self.flowArea(Fi)/6/speed
        return self.dm_record

    def connect_to(self,last_component,next_component):
        self.last=last_component
        last_component.next=self

        self.next=next_component
        next_component.last=self

    def plotValveLift(self, anotherValve=None):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 10))

        tempdata = [mod(each) for each in self.valve_lift.table[0].data]
        ax.plot(self.valve_lift.table[0].data, self.valve_lift.table[1].data)
        ax.set_ylabel(self.valve_lift.table[1].ColName + "(" + self.valve_lift.table[1].ColUnit + ")")
        ax.set_xlabel(self.valve_lift.table[0].ColName + "(" + self.valve_lift.table[0].ColUnit + ")")
        ax.set_xlim([-540, 540])
        ax.set_ylim([0, 1.1 * max(self.valve_lift.table[1].data)])

        # 做标记
        ax.annotate('%.2f $^\circ$ CA' % self.valve_lift.table[0].data[0],
                    xy=(self.valve_lift.table[0].data[0], self.valve_lift.table[1].data[0]), xycoords='data',
                    xytext=(-80, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.annotate('%.2f $^\circ$ CA' % self.valve_lift.table[0].data[-1],
                    xy=(self.valve_lift.table[0].data[-1], self.valve_lift.table[1].data[-1]), xycoords='data',
                    xytext=(0, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        i = self.valve_lift.findMaxValueIndex(1)
        ax.annotate('%.5f' % self.valve_lift.table[1].data[i],
                    xy=(self.valve_lift.table[0].data[i], self.valve_lift.table[1].data[i]), xycoords='data',
                    xytext=(0, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.valve_lift.table[0].data[i], self.valve_lift.table[1].data[i], color="r")
        ax.scatter(self.valve_lift.table[0].data[-1], self.valve_lift.table[1].data[-1], color="r")
        ax.scatter(self.valve_lift.table[0].data[0], self.valve_lift.table[1].data[0], color="r")

        if anotherValve is not None:
            ax.plot(anotherValve.valve_lift.table[0].data, anotherValve.valve_lift.table[1].data)
            ax.set_ylim([0, 1.1 * max(anotherValve.valve_lift.table[1].data)])
            # 做标记
            ax.annotate('%.2f $^\circ$ CA' % anotherValve.valve_lift.table[0].data[0],
                        xy=(anotherValve.valve_lift.table[0].data[0], anotherValve.valve_lift.table[1].data[0]),
                        xycoords='data',
                        xytext=(-80, 30), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.annotate('%.2f $^\circ$ CA' % anotherValve.valve_lift.table[0].data[-1],
                        xy=(anotherValve.valve_lift.table[0].data[-1], anotherValve.valve_lift.table[1].data[-1]),
                        xycoords='data',
                        xytext=(0, 30), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            i = anotherValve.valve_lift.findMaxValueIndex(1)
            ax.annotate('%.5f' % anotherValve.valve_lift.table[1].data[i],
                        xy=(anotherValve.valve_lift.table[0].data[i], anotherValve.valve_lift.table[1].data[i]),
                        xycoords='data',
                        xytext=(0, 30), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.scatter(anotherValve.valve_lift.table[0].data[i], anotherValve.valve_lift.table[1].data[i], color="r")
            ax.scatter(anotherValve.valve_lift.table[0].data[-1], anotherValve.valve_lift.table[1].data[-1], color="r")
            ax.scatter(anotherValve.valve_lift.table[0].data[0], anotherValve.valve_lift.table[1].data[0], color="r")

        # 画垂直线
        ax.axvline(x=0, color='g', linestyle=":")
        ax.axvline(x=180, color='g', linestyle=":")
        ax.axvline(x=-180, color='g', linestyle=":")
        ax.axvline(x=360, color='g', linestyle=":")
        ax.axvline(x=-360, color='g', linestyle=":")

        plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])
        plt.tight_layout()
        plt.show()


class ValveDesign:
    def __init__(self, Bore=130, EVO=55, EVC=15, IVO=20, IVC=40, valvetype="inlet", TDC=0):
        if valvetype == "inlet":
            self.ValveDiameter = 0.39 * Bore
        elif valvetype == "Exhaust":
            self.ValveDiameter = 0.32 * Bore
        else:
            raise Exception("Valve type error, only inlet and exhaust are allowed")

        self.TDC = TDC

        self.EVO = mod(TDC + 180 - EVO, TDC)
        self.EVC = mod(TDC + 360 + EVC - 720, TDC)
        self.IVO = mod(TDC - 360 - IVO, TDC)
        self.IVC = mod(TDC - 180 + IVC, TDC)
        print("EVO={}".format(self.EVO))
        print("EVC={}".format(self.EVC))
        print("IVO={}".format(self.IVO))
        print("IVC={}".format(self.IVC))

    def changeTiming(self, EVO=None, EVC=None, IVO=None, IVC=None, TDC=0):
        if EVO is not None:
            self.EVO = mod(EVO, TDC)
        if EVC is not None:
            self.EVC = mod(EVC, TDC)
        if IVO is not None:
            self.IVO = mod(IVO, TDC)
        if IVC is not None:
            self.IVC = mod(IVC, TDC)

    def changeTDC(self, TDC=360):
        self.EVO += TDC
        self.EVC += TDC
        self.IVO += TDC
        self.IVC += TDC

    def plot(self):
        from ArrayTable import ArrayTable
        table = ArrayTable()
        table.readCSVFile("CylinderPressure.csv")
        for i in range(table.row):
            table.table[0].data[i] = mod(table.table[0].data[i], TDC=self.TDC)
        table.doQuickSort()

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 5))
        ax.plot(table.table[0].data, table.table[1].data)
        plt.xlabel(table.table[0].ColName + "(" + table.table[0].ColUnit + ")")
        plt.ylabel(table.table[1].ColName + "(" + table.table[1].ColUnit + ")")

        ax.scatter(self.IVC, table.linearInterpolate(self.IVC, 1))
        ax.annotate('IVC %.3g $^\circ$CA' % self.IVC,
                    xy=(self.IVC, table.linearInterpolate(self.IVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.IVO, table.linearInterpolate(self.IVO, 1))
        ax.annotate('IVO %.3g $^\circ$CA' % self.IVO,
                    xy=(self.IVO, table.linearInterpolate(self.IVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.EVO, table.linearInterpolate(self.EVO, 1))
        ax.annotate('EVO %.3g $^\circ$CA' % self.EVO,
                    xy=(self.EVO, table.linearInterpolate(self.EVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(self.EVC, table.linearInterpolate(self.EVC, 1))
        ax.annotate('EVC %.3g $^\circ$CA' % self.EVC,
                    xy=(self.EVC, table.linearInterpolate(self.EVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))
        plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])

        # ax.axhline(y=0, color='r', linestyle="-.")
        ax.axvline(x=0, color='g', linestyle=":")
        ax.axvline(x=180, color='g', linestyle=":")
        ax.axvline(x=-180, color='g', linestyle=":")
        ax.axvline(x=360, color='g', linestyle=":")
        ax.axvline(x=-360, color='g', linestyle=":")

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # temp = ValveDesign()
    # temp.plot()

    valvelift = ArrayTable()
    valvelift.readCSVFile("IntakeValveLift.csv")
    val1 = Valve(120e-3, -370, -170, valvelift)
    val1.valve_lift.plot()
    val1.flowAreaData.plot()
    val1.flow_coeff.plot()

    val2 = Valve(12e-3, 150, 370, valvelift, 0, 1.5)
    val2.valve_lift.plot()
    val1.plotValveLift(val2)
