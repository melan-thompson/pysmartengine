#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/9/16
# @Author  : github.com/Melan_Thompson

from .Table import ArrayTable
from .GasProperty import *
from abc import ABCMeta, ABC
from abc import abstractmethod
from math import pi


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
    if abs(TMP2 - 1) < 1.e-5:
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


# 阀门升程曲线抽象类
class ValveLift(metaclass=ABCMeta):
    def __init__(self):
        self.data = ArrayTable(2, 0)
        self.data.setTableUnit(["CA", "m"])
        self.data.setTableHeader(["Crank angle", "Valve lift"])

    @abstractmethod
    def lift(self, crank_angle):
        pass

    @abstractmethod
    def max_lift(self):
        pass

    @abstractmethod
    def plot(self):
        pass

    @abstractmethod
    def move_timing(self, opened, closed):
        pass


# 阀门升程曲线抽象类
class ValveABC(metaclass=ABCMeta):
    def __init__(self):
        # 记录流量的变化
        self.data = ArrayTable(3, 0)
        self.data.setTableHeader(["time", "Pressure ratio", "Mass flow rate"])
        self.data.setTableUnit(["s", "/", "kg/s"])

        # 设置当前时刻
        self.t = None

        # 连接关系
        self.next = None
        self.last = None

        self.dm_record = None

    # 连接关系
    @abstractmethod
    def connect(self, connect_from, connect_to):
        pass

    # 当前时刻的流量
    @abstractmethod
    def mass_flow_rate(self):
        pass

    # 某时刻的有效流通面积
    @abstractmethod
    def flow_area(self):
        pass

    # 更新t时刻的流量值
    @abstractmethod
    def update(self):
        pass


class ValveLiftFunc(ValveLift, ABC):
    def __init__(self, x_max=8.e-3, tau=10):
        super(ValveLiftFunc, self).__init__()
        self.LiftMax = x_max
        self.tau = tau

        def fun(crank_angle):
            from math import cos, pi
            return x_max * (0.5 - 0.5 * cos(2. * crank_angle / tau * pi))

        from numpy import arange
        for i in arange(0, tau, tau / 200):
            self.data.append([i, fun(i)])

    def lift(self, crank_angle):
        if crank_angle<self.data.table[0].data[0] or crank_angle>self.data.table[0].data[-1]:return 0
        else:
            return self.data.linearInterpolate(crank_angle,1)
        # from math import cos, pi
        # return self.LiftMax * (0.5 - 0.5 * cos(2. * (crank_angle-self) / self.tau * pi))

    def max_lift(self):
        return max(self.data.table[1].data)

    def move_timing(self, opened, closed):
        self.tau=closed-opened
        self.data.table[0].data = self.data.table[0].scaleToRange([opened, closed])

    def plot(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 10))
        from numpy import array
        ax.plot(self.data.table[0].data, array(self.data.table[1].data) * 1.e3)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Valve lift(mm)")
        plt.tight_layout()
        plt.show()


class ValveSimple(ValveABC, ABC):
    def __init__(self, _flowArea=1.):
        super(ValveSimple, self).__init__()
        self.flowArea = _flowArea

    # def __massFlowRate(self, p1, T1, R1, k1, p2, T2=0, R2=0, k2=0):
    #     return self.flowArea * flowUnitArea(p1, T1, R1, k1, p2, T2, R2, k2)

    def mass_flow_rate(self):
        return self.flowArea * flowUnitArea(self.last.p, self.last.T, self.last.Rg, self.last.k, self.next.p,
                                            self.next.T, self.next.Rg, self.next.k)

    def connect(self, connect_from, connect_to):
        self.last = connect_from
        self.last.next = self

        self.next = connect_to
        self.next.last = self

    def update(self):
        self.dm_record = self.mass_flow_rate()
        self.record(self.t)
        return self.dm_record

    def record(self, t):
        self.data.append([t, self.next.p / self.last.p, self.dm_record])

    def flow_area(self):
        return self.flowArea

    # def airFlowExample(self, p1, T1, p2=None, T2=None, plot=True):
    #     # from GasProperty import Rg, k_Justi
    #     R1 = Rg()
    #     R2 = Rg()
    #     k1 = k_Justi(T1)
    #     step = 1.e-3 * p1
    #     import numpy as np
    #     table = ArrayTable(2, 0)
    #     for i in np.arange(0, p1 + step, step):
    #         table.append([i / p1, self.__massFlowRate(p1, T1, R1, k1, i)])
    #
    #     if p2 is not None and T2 is not None:
    #         k2 = k_Justi(T2)
    #         step2 = 0.001 * p2
    #         for i in np.arange(0.4 * p2, p2, step2):
    #             table.append([p2 / i, self.__massFlowRate(i, T1, R1, k1, p2, T2, R2, k2)])
    #
    #     table.doQuickSort()
    #     if plot is True:
    #         import matplotlib.pyplot as plt
    #         plt.plot(table.table[0].data, table.table[1].data)
    #         choke1 = pow(2 / (k1 + 1), k1 / (k1 - 1))
    #         print(choke1)
    #         massflow1 = self.__massFlowRate(p1, T1, R1, k1, p1 * choke1)
    #         plt.scatter(choke1, massflow1)
    #
    #         choke2 = pow(2 / (k2 + 1), k2 / (k2 - 1))
    #         print(1 / choke2)
    #         massflow2 = self.__massFlowRate(p1 * choke2, T1, R1, k1, p2, T2, R2, k2)
    #         plt.scatter(1 / choke2, massflow2)
    #
    #         plt.axvline(x=1, color='r', linestyle=":")
    #         plt.axhline(y=0, color='r', linestyle=":")
    #
    #         plt.annotate('%.5f' % choke1,
    #                      xy=(choke1, massflow1), xycoords='data',
    #                      xytext=(-80, -30), textcoords='offset points',
    #                      arrowprops=dict(arrowstyle="->"))
    #
    #         plt.annotate('%.5f' % (1 / choke2),
    #                      xy=(1 / choke2, massflow2), xycoords='data',
    #                      xytext=(40, 30), textcoords='offset points',
    #                      arrowprops=dict(arrowstyle="->"))
    #         plt.xlabel("Pressure ratio(/)")
    #         plt.ylabel("Mass flow rate(kg/s)")
    #
    #         plt.tight_layout()
    #         plt.show()
    #
    #     return table


class Valve(ValveABC, ABC):
    def __init__(self, valveDiameter, valve_lift, sigma=0, _valve_flow_coeff=None):
        super(Valve, self).__init__()

        self.valve_lift = valve_lift
        self.valve_diameter = valveDiameter
        self.sigma = sigma * pi / 180.

        # 检查阀门的直径是否正确
        maxlift = self.valve_lift.max_lift()
        print("maximum valve lift is {}".format(maxlift))
        if valveDiameter < maxlift / 0.544:
            raise Exception("valve diameter is too small, it can not be less than {}".format(maxlift / 0.544))

        # 流量系数
        self.flow_coeff = _valve_flow_coeff
        if _valve_flow_coeff is None:
            self.flow_coeff = ArrayTable(2, 0)
            self.flow_coeff.setTableHeader(["Valve lift", "Flow coefficient"])
            self.flow_coeff.setTableUnit(["m", "/"])
            import numpy as np
            for i in np.arange(0, maxlift, maxlift / 100):
                flowcoeff = self.__flowcoefficient(i)
                if flowcoeff < 0:
                    raise Exception("valve flow coefficient can not less than 0")
                self.flow_coeff.append([i, flowcoeff])

        # 有效流通面积，先计算好，不用再算
        self.flowAreaData = ArrayTable(2, 0)
        self.flowAreaData.setTableHeader(["Valve lift", "Flow area"])
        self.flowAreaData.setTableUnit(["$^\circ$CA", "$m^2$"])
        from numpy import arange
        for i in arange(-540, 540):
            self.t = i
            self.flowAreaData.append([i, self.flow_area()])
        self.t = None

    # 推荐的流量系数
    def __flowcoefficient(self, _valvelift):
        return 0.98 - 3.3 * pow((_valvelift / self.valve_diameter), 2)

    # 计算当前时刻的有效流通面积
    def flow_area(self):
        from math import pi, cos, sin
        hv = self.valve_lift.lift(self.t)
        flowcoeff = self.__flowcoefficient(hv)
        if flowcoeff < 0:
            raise Exception("valve flow coefficient can not less than 0")
        return flowcoeff * pi * hv * cos(self.sigma) * (self.valve_diameter + hv * sin(self.sigma) * cos(self.sigma))

    def connect(self, connect_from, connect_to):
        self.last = connect_from
        connect_from.next = self

        self.next = connect_to
        connect_to.last = self

    # 计算质量流量
    def mass_flow_rate(self):
        self.dm_record = flowUnitArea(self.last.p, self.last.T, self.last.Rg, self.last.k, self.next.p,
                                      self.next.T, self.next.Rg, self.next.k) * self.flow_area()
        return self.dm_record

    def update(self):
        pass

    # def plotValveLift(self, anotherValve=None):
    #     import matplotlib.pyplot as plt
    #     fig, ax = plt.subplots(1, figsize=(10, 10))
    #
    #     tempdata = [mod(each) for each in self.valve_lift.table[0].data]
    #     ax.plot(self.valve_lift.table[0].data, self.valve_lift.table[1].data)
    #     ax.set_ylabel(self.valve_lift.table[1].ColName + "(" + self.valve_lift.table[1].ColUnit + ")")
    #     ax.set_xlabel(self.valve_lift.table[0].ColName + "(" + self.valve_lift.table[0].ColUnit + ")")
    #     ax.set_xlim([-540, 540])
    #     ax.set_ylim([0, 1.1 * max(self.valve_lift.table[1].data)])
    #
    #     # 做标记
    #     ax.annotate('%.2f $^\circ$ CA' % self.valve_lift.table[0].data[0],
    #                 xy=(self.valve_lift.table[0].data[0], self.valve_lift.table[1].data[0]), xycoords='data',
    #                 xytext=(-80, 30), textcoords='offset points',
    #                 arrowprops=dict(arrowstyle="->"))
    #
    #     ax.annotate('%.2f $^\circ$ CA' % self.valve_lift.table[0].data[-1],
    #                 xy=(self.valve_lift.table[0].data[-1], self.valve_lift.table[1].data[-1]), xycoords='data',
    #                 xytext=(0, 30), textcoords='offset points',
    #                 arrowprops=dict(arrowstyle="->"))
    #
    #     i = self.valve_lift.findMaxValueIndex(1)
    #     ax.annotate('%.5f' % self.valve_lift.table[1].data[i],
    #                 xy=(self.valve_lift.table[0].data[i], self.valve_lift.table[1].data[i]), xycoords='data',
    #                 xytext=(0, 30), textcoords='offset points',
    #                 arrowprops=dict(arrowstyle="->"))
    #
    #     ax.scatter(self.valve_lift.table[0].data[i], self.valve_lift.table[1].data[i], color="r")
    #     ax.scatter(self.valve_lift.table[0].data[-1], self.valve_lift.table[1].data[-1], color="r")
    #     ax.scatter(self.valve_lift.table[0].data[0], self.valve_lift.table[1].data[0], color="r")
    #
    #     if anotherValve is not None:
    #         ax.plot(anotherValve.valve_lift.table[0].data, anotherValve.valve_lift.table[1].data)
    #         ax.set_ylim([0, 1.1 * max(anotherValve.valve_lift.table[1].data)])
    #         # 做标记
    #         ax.annotate('%.2f $^\circ$ CA' % anotherValve.valve_lift.table[0].data[0],
    #                     xy=(anotherValve.valve_lift.table[0].data[0], anotherValve.valve_lift.table[1].data[0]),
    #                     xycoords='data',
    #                     xytext=(-80, 30), textcoords='offset points',
    #                     arrowprops=dict(arrowstyle="->"))
    #
    #         ax.annotate('%.2f $^\circ$ CA' % anotherValve.valve_lift.table[0].data[-1],
    #                     xy=(anotherValve.valve_lift.table[0].data[-1], anotherValve.valve_lift.table[1].data[-1]),
    #                     xycoords='data',
    #                     xytext=(0, 30), textcoords='offset points',
    #                     arrowprops=dict(arrowstyle="->"))
    #
    #         i = anotherValve.valve_lift.findMaxValueIndex(1)
    #         ax.annotate('%.5f' % anotherValve.valve_lift.table[1].data[i],
    #                     xy=(anotherValve.valve_lift.table[0].data[i], anotherValve.valve_lift.table[1].data[i]),
    #                     xycoords='data',
    #                     xytext=(0, 30), textcoords='offset points',
    #                     arrowprops=dict(arrowstyle="->"))
    #
    #         ax.scatter(anotherValve.valve_lift.table[0].data[i], anotherValve.valve_lift.table[1].data[i], color="r")
    #         ax.scatter(anotherValve.valve_lift.table[0].data[-1], anotherValve.valve_lift.table[1].data[-1], color="r")
    #         ax.scatter(anotherValve.valve_lift.table[0].data[0], anotherValve.valve_lift.table[1].data[0], color="r")
    #
    #     # 画垂直线
    #     ax.axvline(x=0, color='g', linestyle=":")
    #     ax.axvline(x=180, color='g', linestyle=":")
    #     ax.axvline(x=-180, color='g', linestyle=":")
    #     ax.axvline(x=360, color='g', linestyle=":")
    #     ax.axvline(x=-360, color='g', linestyle=":")
    #
    #     plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])
    #     plt.tight_layout()
    #     plt.show()


class ValveDesign:
    """
    气门正时类
    """
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
    V = ValveLiftFunc(14e-3, tau=190)
    V.move_timing(150, 300)
    V.plot()

    # temp = ValveDesign()
    # temp.plot()

    v = ValveSimple()
    v.airFlowExample(1.e5, 300, 1.5e5, 400, plot=True)

    valvelift = ArrayTable()
    valvelift.readCSVFile("IntakeValveLift.csv")
    val1 = Valve(120e-3, -370, -170, valvelift)
    val1.valve_lift.plot()
    val1.flowAreaData.plot()
    val1.flow_coeff.plot()

    val2 = Valve(12e-3, 150, 370, valvelift, 0, 1.5)
    val2.valve_lift.plot()
    val1.plotValveLift(val2)
