import sys

sys.path.append("../Properties of working fluids/")
from GasProperty import *

def twoPointInterpolate(x, _x=[0, 1], _y=[2, 3]):
    def fun(x):
        return (x - _x[0]) / (_x[1] - _x[0]) * _y[1] + (x - _x[1]) / (_x[0] - _x[1]) * _y[0]

    return fun(x)


def TurbochargePressure(pme=10e5, T_im=300, eta_et=0.36, phi_a=1.1, VE=0.86, L0=14.3, Hu=42700e3):
    """
    计算增压压力
    :param pme: 平均有效压力(Pa)
    :param T_im: 进气管温度(K)
    :param eta_et: 有效热效率
    :param phi_a: 过量空气系数
    :param VE: 充量系数，自然吸气式发动机最大为90%
    :param L0: 燃空当量比
    :param Hu: 燃料低热值
    :return: 增压压力(Pa)
    """
    from GasProperty import Rg
    return pme * T_im / eta_et * phi_a / VE * Rg() * L0 / Hu


def TAfterCompressor(T0, pik, etak=1, tau=1):
    """
    计算压气机后的温度
    :param T0: 环境温度(K)
    :param pik: 压气机压比
    :param etak: 压气机效率,等于1时为等熵压缩温度
    :param tau: 考虑向外散热的冷却系数，tau=1.04~1.10
    :return:压气机后的温度(K)
    """
    k = k_Justi(T0)
    result = T0 + T0 * (pow(pik, (k - 1) / k) - 1) / (etak * tau)
    return result


def TAfterTurbine(Tt, pit, etat=1):
    k = k_exhaust_gu(Tt)
    return Tt * (1 - etat * (1 - pow(pit, -(k - 1) / k)))


def piK(Tt, pit, etat, etak, T0=273.15 + 20, air_fuel_ratio=23.85, ma=1, me=None, etam=1):
    """
    由涡轮前的状态和质量流量计算压气机的压比
    :param Tt:压气机前的温度(K)
    :param pit:涡轮压降
    :param etat:涡轮效率
    :param etak:压气机效率
    :param air_fuel_ratio:空燃比
    :param T0:环境温度(K)
    :param ma:空气质量流量，计算压气机功率时用到的，否则不用
    :param me:废气质量流量，计算涡轮功率用到的，否则不用
    :param etam:涡轮增压器机械效率
    :return:压气机压比
    """
    if me is None:me=ma+ma/air_fuel_ratio

    gammak = k_Justi(T0)
    gammat = k_exhaust_gu(Tt)
    print("Turbine power:{}kW".format(
        1.e-3 * me * cp_exhaust_gu(Tt) * Tt * (1 - pow(pit, -(gammat - 1) / gammat)) * etat))
    print("Temperature after turbine:{}K".format(TAfterTurbine(Tt, pit, etat)))

    def fun(pik):
        return air_fuel_ratio * cp_Justi(T0) * T0 * (pow(pik, (gammak - 1) / gammak) - 1) - (
                1 + air_fuel_ratio) * cp_exhaust_gu(Tt) * Tt * (
                       1 - pow(pit, -(gammat - 1) / gammat)) * etat * etak * etam

    result = 2
    h = 1.e-3
    while abs(fun(result)) > 1.e-7:
        result -= fun(result) / ((fun(result + h) - fun(result - h)) / (2 * h))

    print("Compressor power:{}kW".format(
        1.e-3 * ma * cp_Justi(T0) * T0 * (pow(result, (gammak - 1) / gammak) - 1) / etak))
    print("Temperature after compressor:{}K".format(TAfterCompressor(T0, result, etak)))

    return result


def piT(pik, etatk=1, Tt=600, T0=273.15 + 20, air_fuel_ratio=23.85):
    gammak = k_Justi(T0)
    gammat = k_exhaust_gu(Tt)
    temp = (cp_Justi(T0) * T0 * air_fuel_ratio) / ((air_fuel_ratio + 1) * cp_exhaust_gu(Tt) * Tt * etatk)
    temp2 = pow(pik, (gammak - 1) / gammak) - 1
    return pow(1 - temp * temp2, -gammat / (gammat - 1))


class Compressor:
    def __init__(self, mapfile, Tref=298, Pref=101325):
        from pandas import read_excel
        self.map = read_excel(mapfile)
        self.Tref = Tref
        self.Pref = Pref

        # 获取表头
        self.headers = list(self.map.columns)

        # 对转速、流量进行排序
        self.map.sort_values(by=[self.headers[0], self.headers[1]])

        # 按照转速分组
        self.group = self.map.groupby(self.headers[0])

        # 获取转速
        self.speeds = list()
        for each, ii in self.group:
            self.speeds.append(each)

    def plot(self, showlegend=False, hightlight=None):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        for each in self.speeds:
            ax[1] = self.group.get_group(each).plot(x=self.headers[1], y=self.headers[2],
                                                    label="speed=" + str(each), ax=ax[1], marker="x")
            ax[0] = self.group.get_group(each).plot(x=self.headers[1], y=self.headers[3], label="speed=" + str(each),
                                                    ax=ax[0], marker="x")
        # surge line
        xx = [];yy = []
        for each in self.speeds:
            xx.append(self.group.get_group(each).iloc[0, 1])
            yy.append(self.group.get_group(each).iloc[0, 2])
        ax[1].plot(xx, yy, "r--")

        # choke flow line
        xx1 = [];
        yy1 = []
        for each in self.speeds:
            xx1.append(self.group.get_group(each).iloc[-1, 1])
            yy1.append(self.group.get_group(each).iloc[-1, 2])
        ax[1].plot(xx1, yy1, "r--")

        # 最大效率线
        # xx2=[]
        # yy2=[]
        # zz2=[]
        # for each in self.speeds:
        #     table=self.group.get_group(each)
        #     indexx=table[table[self.headers[3]] == max(table[self.headers[3]])].index[0]
        #     print(indexx)
        #     xx2.append(self.map.iloc[indexx,1])
        #     yy2.append(self.map.iloc[indexx,2])
        #     zz2.append(self.map.iloc[indexx,3])
        # ax[0].plot(xx2,zz2,'r--')
        # ax[1].plot(xx2,yy2,'r--')

        ax[0].set_ylabel(self.headers[3])
        ax[1].set_ylabel(self.headers[2])

        ax[1].annotate('surge line',
                       xy=(xx[(len(xx) // 2)], yy[(len(xx) // 2)]), xycoords='data',
                       xytext=(-50, 30), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->"))
        ax[1].annotate('choke flow line',
                       xy=(xx1[(len(xx) // 2)], yy1[(len(xx) // 2)]), xycoords='data',
                       xytext=(-10, -40), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->"))

        # 删除图例
        if showlegend:
            plt.legend(loc='best')
        else:
            ax[0].get_legend().remove()
            ax[1].get_legend().remove()

        if hightlight is not None:
            ax[1].scatter(hightlight[0], hightlight[1], marker='*', c='r')

        plt.tight_layout()
        plt.show()

    def addLine(self, speed):
        if speed < min(self.speeds) or speed > max(self.speeds):
            raise Exception("interpolate out of range")
        import pandas as pd
        newgroup = pd.DataFrame(index=self.group.get_group(self.speeds[0]).index,
                                columns=self.group.get_group(self.speeds[0]).columns)
        newgroup["corrected speed,RPM"] = speed

        for i in range(len(self.speeds)):
            if self.speeds[i] >= speed:
                speedleft = self.speeds[i - 1]
                speedright = self.speeds[i]
                break

        for j in range(len(self.group.get_group(self.speeds[0]))):
            for k in range(1, 4):
                temp1 = self.group.get_group(speedleft).iloc[j, k]
                temp2 = self.group.get_group(speedright).iloc[j, k]
                newgroup.loc[j, self.headers[k]] = twoPointInterpolate(speed, [speedleft, speedright], [temp1, temp2])

        self.map = self.map.append(newgroup)

        # 对转速排序
        self.map.sort_values(by=[self.headers[0],self.headers[1]])

        # 分组
        self.group = self.map.groupby(self.headers[0])

        # 获取转速
        self.speeds = list()
        for each, ii in self.group:
            self.speeds.append(each)

    def interpolate(self, order=3):
        from numpy import linspace
        speeds = self.speeds
        for i in range(len(self.speeds) - 1):
            space = linspace(speeds[i], speeds[i + 1], order + 2)[1:-1]
            for each in space:
                self.addLine(each)

    def massFlowrate(self, speed, pik):
        if speed not in self.speeds:
            self.addLine(speed)

        # 压比
        temp = self.group.get_group(speed)[self.headers[2]]
        # 流量
        temp2 = self.group.get_group(speed)[self.headers[1]]

        # 效率
        temp3 = self.group.get_group(speed)[self.headers[3]]

        for i in range(len(temp)):
            if pik < min(temp):
                raise Exception("pressure ratio is too small while calculating mass flow rate")
            elif pik > max(temp):
                raise Exception("pressure ratio is too big while calculating mass flow rate")
            if temp[i] < pik:
                massflowrate = twoPointInterpolate(pik, [temp[i - 1], temp[i]], [temp2[i - 1], temp2[i]])
                break

        for j in range(len(temp)):
            if massflowrate < temp2[j]:
                print(j)
                efficiency = twoPointInterpolate(massflowrate, [temp2[j - 1], temp2[j]], [temp3[j - 1], temp3[j]])
                break
        from math import sqrt
        print("Speed={},pressure ratio={},mass flow rate={}kg/s,efficiency={}".format(speed*sqrt(self.Tref), pik, massflowrate*self.Pref/sqrt(self.Tref), efficiency))
        return massflowrate


class Turbine:
    def __init__(self, mapfile):
        from pandas import read_excel
        self.map = read_excel(mapfile)
        # self.Tref = Tref
        # self.Pref = Pref

        # 获取表头
        self.headers = list(self.map.columns)

        # 对转速、流量进行排序
        self.map.sort_values(by=[self.headers[0], self.headers[1]])

        # 按照转速分组
        self.group = self.map.groupby(self.headers[0])

        # 获取转速
        self.speeds = list()
        for each, ii in self.group:
            self.speeds.append(each)

    def plot(self, showlegend=False, hightlight=None):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        for each in self.speeds:
            ax[1] = self.group.get_group(each).plot(x=self.headers[2], y=self.headers[1],
                                                    label="speed=" + str(each), ax=ax[1], marker="x")
            ax[0] = self.group.get_group(each).plot(x=self.headers[1], y=self.headers[3], label="speed=" + str(each),
                                                    ax=ax[0], marker="x")

        ax[0].set_ylabel(self.headers[3])
        ax[1].set_ylabel(self.headers[2])

        # 删除图例
        if showlegend:
            plt.legend(loc='best')
        else:
            ax[0].get_legend().remove()
            ax[1].get_legend().remove()

        if hightlight is not None:
            ax[1].scatter(hightlight[0], hightlight[1], marker='*', c='r')

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    C = Compressor("./Data/CompressorMapMa.xlsx")
    C.plot()
    C.interpolate(3)
    # C.massFlowrate(60000,0)
    C.plot(showlegend=True)
    hightlight = [C.massFlowrate(75000, 2.4), 2.4]

    Turb = Turbine("./Data/TurbineMapMa.xlsx")
    Turb.plot(showlegend=True)
