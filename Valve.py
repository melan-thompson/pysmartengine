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
    from math import *
    TMP2 = p2 / p1
    if TMP2 < 1.0:  ##正向流动
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


class ValveSimple:
    def __init__(self, _flowArea=1.):
        from ArrayTable import ArrayTable
        self.flowArea = 1
        self.data = ArrayTable(2, 0)
        self.data.setTableHeader(["Pressure ratio", "Mass flow rate"])
        self.data.setTableUnit(["/", "kg/s"])

    def __massFlowRate(self, p1, T1, R1, k1, p2, T2=0, R2=0, k2=0):
        return self.flowArea * flowUnitArea(p1, T1, R1, k1, p2, T2, R2, k2)

    def airFlowExample(self, p1, T1, p2=None, T2=None, K2=None):
        from GasProperty import *
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
        open = _valve_lift.table[0].data[0]
        close = _valve_lift.table[0].data[-1]
        temp = (close_timing_angle - open_timing_angle) / (close - open)
        self.valve_lift = ArrayTable(2, 0)
        self.valve_lift.setTableUnit(["CA", "m"])
        self.valve_lift.setTableHeader(["Crank angle", "Valve lift"])
        for i in range(_valve_lift.row):
            self.valve_lift.append(
                [self.open_timing_angle + temp * (_valve_lift.table[0].data[i] - _valve_lift.table[0].data[0]),
                 _valve_lift.table[1].data[i] * rock])

        self.flow_coeff = _valve_flow_coeff
        if _valve_flow_coeff is None:
            self.flow_coeff = ArrayTable(2, 0)
            self.flow_coeff.setTableHeader(["Crank angle", "Flow coefficient"])
            self.flow_coeff.setTableUnit(["CA", "/"])
            maxlift = self.valve_lift.findMaxValue(1)
            import numpy as np
            for i in np.arange(0, maxlift, 0.001):
                self.flow_coeff.append([i, self.__flowcoefficient(i)])

    def valveLiftStartFromZero(self):
        temp = self.valve_lift
        startpoint = self.valve_lift.table[0].data[0]
        temp.table[0] + (-startpoint)
        return temp

    def __flowcoefficient(self, _valvelift):
        return 0.98 - 3.3 * pow((_valvelift / self.valve_diameter), 2)
