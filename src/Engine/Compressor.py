#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/9/16
# @Author  : github.com/Melan_Thompson

from .GasProperty import *
from .Table import *
from abc import ABCMeta, ABC
from abc import abstractmethod

# 两点插值
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
    from .GasProperty import k_Justi
    k = k_Justi(T0)
    result = T0 + T0 * (pow(pik, (k - 1) / k) - 1) / (etak * tau)
    return result


def TAfterTurbine(Tt, pit, etat=1):
    """
    计算涡轮后温度
    :param Tt: 涡轮前温度,K
    :param pit: 涡轮膨胀比
    :param etat: 涡轮等熵效率
    :return: 涡轮后温度
    """
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
    if me is None: me = ma + ma / air_fuel_ratio

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


def ETpiK(power, etak=0.5, Tbefore=273.15):
    pass


def piT(pik, etatk=1, Tt=600, T0=273.15 + 20, air_fuel_ratio=23.85):
    # \phi\frac{{{c_{pa}}{T_0}\left( {{\pi _k}^{\frac{{\gamma  - 1}}{\gamma }} - 1} \right)}}{{{\eta _k}}} = \left( {\phi  + 1} \right){c_{pe}}{T_t}\left( {1 - {\pi _t}^{ - \frac{{\gamma  - 1}}{\gamma }}} \right){\eta _t}{\eta _m}
    """
    计算涡轮膨胀比
    :param pik: 压气机压比
    :param etatk: 涡轮增压器总效率
    :param Tt: 涡轮前温度,K
    :param T0: 压气机前温度,K
    :param air_fuel_ratio: 空燃比
    :return: 涡轮膨胀比
    """
    gammak = k_Justi(T0)
    gammat = k_exhaust_gu(Tt)
    temp = (cp_Justi(T0) * T0 * air_fuel_ratio) / ((air_fuel_ratio + 1) * cp_exhaust_gu(Tt) * Tt * etatk)
    temp2 = pow(pik, (gammak - 1) / gammak) - 1
    return pow(1 - temp * temp2, -gammat / (gammat - 1))


class CompressorSimple():
    def __init__(self,T0=273.15+30,p0=1.e5) -> None:
        self.T0=T0
        self.p0=p0
        self.data=ArrayTable(4,0)
        self.data.setTableHeader(["Compression ratio","Mass flow rate","efficiency","power"])
        self.data.setTableUnit(["/","kg/s","/","W"])

        self.dm=None

        # 连接关系
        self.last=None
        self.next=None

    # 计算压气机后的温度
    def Tk(self,eff,pik):
        # 计算压气机后的温度
        return TAfterCompressor(self.T0,pik,eff)

    # 由流量计算压气机功
    def power(self,mass):
        pass

    def _mass(self,power,pik,eff):
        """
        由功计算增压器流量
        """
        Ttemp=self.Tk(eff,pik)
        self.dm=power/(cp_Justi(Ttemp)*Ttemp-cp_Justi(self.T0)*self.T0)
        return self.dm

    def connect_to(self,connect_from,connect_to):
        connect_from.next=self
        connect_from.last=self

        self.next=connect_to
        self.last=connect_from

    def mass(self):
        return 
    # 动态特性和静态特性 合作
    # 气路传递和能量传递,能量导致温度等变化。
        
# 气路port和能量port
class port:
    def __init__(self) -> None:
        pass

    def portNum(self):
        pass


# 涡轮增压器抽象类
class TurbochargerBase(metaclass=ABCMeta):
    def __init__(self,mapfile,Tref,Pref):
        from pandas import read_excel, to_numeric, DataFrame
        self.map = read_excel(mapfile).apply(to_numeric)

        self.Tref = Tref
        self.Pref = Pref

    # 由转速和压比求流量和效率
    @abstractmethod
    def masseff(self,speed,pik):
        pass
    
    # 由转速和压比求实际功率
    @abstractmethod
    def power(self,speed,pik):
        pass
        # mass,eff=self.masseff(speed,pik)

    # 画图
    @abstractmethod
    def plot(self):
        pass




class Compressor:
    def __init__(self, mapfile, Tref=298, Pref=101325, MapType="real", order=2):
        from pandas import read_excel, to_numeric, DataFrame
        self.map = read_excel(mapfile).apply(to_numeric)
        self.Tref = Tref
        self.Pref = Pref
        self.intefunorder = order

        # 等效率表格
        self.effmap = DataFrame(columns=self.map.columns)

        # 获取表头
        self.headers = list(self.map.columns)

        # 对转速、流量进行排序
        self.map = self.map.sort_values(by=[self.headers[0], self.headers[1]])

        # 按照转速分组
        self.group = self.map.groupby(self.headers[0])

        # 获取转速
        self.speeds = list()
        for each, ii in self.group:
            self.speeds.append(each)

        # 加入插值函数{n:[fmass,feff]}
        import scipy.interpolate as inte
        self.intefuns = {}
        for each in self.speeds:
            fmass = inte.interp1d(self.group.get_group(each)[self.headers[1]],
                                  self.group.get_group(each)[self.headers[2]], kind=order)
            feff = inte.interp1d(self.group.get_group(each)[self.headers[1]],
                                 self.group.get_group(each)[self.headers[3]], kind=order)
            funs = {"fmass": fmass, "feff": feff}
            self.intefuns[each] = funs

    def plot(self, showlegend=True, hightlight=None, order=2):
        import matplotlib.pyplot as plt
        from numpy import linspace
        fig, ax = plt.subplots(2, 1, figsize=(20, 10))
        # 添加坐标
        ax[0].set_ylabel(self.headers[3])
        ax[1].set_ylabel(self.headers[2])
        ax[1].set_xlim(0,1.2*max(self.map[self.headers[1]]))

        import scipy.interpolate as inte
        if order >= 2:
            for each in self.speeds:
                xx = linspace(min(self.group.get_group(each)[self.headers[1]]),
                              max(self.group.get_group(each)[self.headers[1]]))
                f1 = inte.interp1d(self.group.get_group(each)[self.headers[1]],
                                   self.group.get_group(each)[self.headers[2]], kind=order)
                f2 = inte.interp1d(self.group.get_group(each)[self.headers[1]],
                                   self.group.get_group(each)[self.headers[3]], kind=order)
                ax[0].plot(xx, f2(xx), label="speed=" + str(round(each,2)))
                ax[1].plot(xx, f1(xx), label="speed=" + str(round(each,2)))

        else:
            for each in self.speeds:
                ax[1] = self.group.get_group(each).plot(x=self.headers[1], y=self.headers[2],
                                                        label="speed=" + str(round(each,2)), ax=ax[1], marker="x")
                ax[0] = self.group.get_group(each).plot(x=self.headers[1], y=self.headers[3],
                                                        label="speed=" + str(round(each,2)),
                                                        ax=ax[0], marker="x")

        # surge line,喘振线
        xx =[]
        yy = []
        for each in self.speeds:
            xx.append(self.group.get_group(each).iloc[0, 1])
            yy.append(self.group.get_group(each).iloc[0, 2])
        maxspeedline=self.group.get_group(max(self.speeds))
        xx1=[i for i in xx]
        yy1=[i for i in yy]
        for i in range(1,len(maxspeedline)):
            xx1.append(maxspeedline.iloc[i, 1])
            yy1.append(maxspeedline.iloc[i, 2])
        self.fsurgeline = inte.interp1d(xx1, yy1, kind=1)
        xx2 = linspace(min(xx), max(xx), 100)
        ax[1].plot(xx2, self.fsurgeline(xx2), "r--",linewidth=4)
        ax[1].annotate('surge line',
                       xy=(xx[(len(xx) // 2)], yy[(len(xx) // 2)]), xycoords='data',
                       xytext=(-50, 30), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->"))

        # choke flow line,堵塞线
        xx =[]
        yy = []
        for each in self.speeds:
            xx.append(self.group.get_group(each).iloc[-1, 1])
            yy.append(self.group.get_group(each).iloc[-1, 2])
        self.fchoke = inte.interp1d(xx, yy, kind=1)
        xx3 = linspace(min(xx), max(xx), 100)
        ax[1].plot(xx3, self.fchoke(xx3), "r--",linewidth=4)
        ax[1].annotate('choke flow line',
                       xy=(xx[(len(xx) // 2)], yy[(len(xx) // 2)]), xycoords='data',
                       xytext=(-10, -40), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->"))

        # # 最大效率线,流量，压比，效率
        xx2 =[]
        yy2 =[]
        zz2 = []
        for each in self.speeds:
            table = self.group.get_group(each)
            indexx = table.idxmax()[3]
            xx2.append(self.group.get_group(each).loc[indexx][1])
            yy2.append(self.group.get_group(each).loc[indexx][2])
            zz2.append(self.group.get_group(each).loc[indexx][3])
        ax[0].plot(xx2, zz2, 'r-.',linewidth=4)
        ax[1].plot(xx2, yy2, 'r-.',linewidth=4)
        ax[1].annotate('maximum efficiency line',
                       xy=(max(xx2), max(yy2)), xycoords='data',
                       xytext=(10, 40), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->"))

        rows=30
        col=30
        # GRP等效率线插值
        x=linspace(min(self.map[self.headers[1]]),max(self.map[self.headers[1]]),col)
        from numpy import array,append,zeros
        lines=linspace(min(self.map[self.headers[3]]),max(self.map[self.headers[3]]),10)
        # print(array([x]*10).ravel())
        xx=array([x]*rows)
        print(xx)
        yy=array([])
        for each in x:
            yy=append(yy,linspace(0,self.fsurgeline(each),rows))
        yy.resize(col,rows)
        yy=yy.transpose()
        print(xx.shape)
        print(yy.shape)

        # 画效率线
        T=ArrayTable(4,0)
        # T.readCSVFile("CompressorMap.csv")
        T.fromPandas(self.map)
        T.GPR([xx[1][1],yy[1][1]],[1,2],3, trained=False)
        z=zeros([rows,col])
        for i in range(rows):
            for j in range(col):
                z[i][j] = T.GPR([xx[i][j], yy[i][j]], [1, 2], 3, trained=True)
        contour=ax[1].contour(xx,yy,z,lines,colors='k')
        ax[1].clabel(contour,fontsize=10,colors=('k','r'))

        # 等效率线
        # print(self.effmap)
        # x = self.effmap[self.headers[1]]
        # y = self.effmap[self.headers[2]]
        # import numpy as np
        # x = np.r_[x, x[0]]
        # y = np.r_[y, y[0]]
        # from scipy import interpolate
        # tck, u = interpolate.splprep([x, y], s=0, per=True)
        # xi, yi = interpolate.splev(np.linspace(0, 1, 1000), tck)
        # ax[1].plot(x, y)

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
                                columns=self.group.get_group(self.speeds[0]).columns, dtype=float)
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

        # 增加插值函数
        import scipy.interpolate as inte
        fmass = inte.interp1d(newgroup[self.headers[1]], newgroup[self.headers[2]],
                              kind=self.intefunorder)
        feff = inte.interp1d(newgroup[self.headers[1]], newgroup[self.headers[3]],
                             kind=self.intefunorder)
        funs = {"fmass": fmass, "feff": feff}
        self.intefuns[speed] = funs

        self.map = self.map.append(newgroup)

        # 对转速排序
        self.map = self.map.sort_values(by=[self.headers[0], self.headers[1]])

        from pandas import to_numeric
        # 重置索引,转化为数值
        self.map = self.map.reset_index(drop=True)

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
        print("maximum pressure ratio={},minimum pressure ratio={} at speed {}".format(max(temp), min(temp), speed))

        # 流量
        temp2 = self.group.get_group(speed)[self.headers[1]]

        # 效率
        temp3 = self.group.get_group(speed)[self.headers[3]]

        if pik < min(temp):
            raise Exception("pressure ratio is too small while calculating mass flow rate")
        elif pik > max(temp):
            raise Exception("pressure ratio is too big while calculating mass flow rate")

        for i in range(len(temp)):
            if temp[i] == pik:
                massflowrate = temp2[i]
                # efficiency=temp3[i]
                break
            elif temp[i] < pik:
                massflowrate = twoPointInterpolate(pik, [temp[i - 1], temp[i]], [temp2[i - 1], temp2[i]])
                break

        for j in range(len(temp)):
            if massflowrate == temp2[j]:
                efficiency = temp3[j]
                break
            if massflowrate < temp2[j]:
                # print(j)
                efficiency = twoPointInterpolate(massflowrate, [temp2[j - 1], temp2[j]], [temp3[j - 1], temp3[j]])
                break
        # from math import sqrt
        # print("Speed={},pressure ratio={},mass flow rate={}kg/s,efficiency={}".format(speed * sqrt(self.Tref), pik,
        #                                                                               massflowrate * self.Pref / sqrt(
        #                                                                                   self.Tref), efficiency))

        return massflowrate, efficiency

    # def pik(self,power,speed):

    def effline(self, eff):
        # 分三种情况，当eff大于最大效率和小于最小效率，直接跳过
        # 只有一个交点，有两个交点
        def listTodic(list1, list2):
            result = {}
            for i in range(len(list1)):
                result[list1[i]] = list2[i]
            return result

        def newto(fun, x0, h=0.01, tol=1.e-5):
            result = x0
            while abs(fun(result)) > tol:
                dfun = (fun(result + h) - fun(result - h)) / 2. / h
                result -= fun(result) / dfun
            return result

        from pandas import DataFrame, merge
        from scipy.optimize import newton
        newgroup = DataFrame(columns=self.map.columns, dtype=float)
        for each in self.speeds:
            if eff > max(self.group.get_group(each)[self.map.columns[3]]) or eff < min(
                    self.group.get_group(each)[self.map.columns[3]]):
                pass
            elif eff < self.group.get_group(each).iloc[0, 3]:
                right = -1
                while self.group.get_group(each).iloc[right, 3] < eff: right -= 1
                y1 = self.group.get_group(each).iloc[right + 1, 3] - eff
                x1 = self.group.get_group(each).iloc[right + 1, 1]
                y2 = self.group.get_group(each).iloc[right, 3] - eff
                x2 = self.group.get_group(each).iloc[right, 1]
                mass = x1 - y1 * (x2 - x1) / (y2 - y1)

                # mass=newton(lambda x:self.intefuns[each]["feff"](x)-eff,self.group.get_group(each)[self.map.columns[1]][-1])
                pi = self.intefuns[each]["fmass"](mass)
                newgroup.append(listTodic(self.map.columns, [each, mass, pi, eff]))
            else:
                # print(self.intefuns[each]["feff"]((self.group.get_group(each).iloc[-1,1]+max(self.group.get_group(each)[self.map.columns[1]]))/2.))
                right = -1
                while self.group.get_group(each).iloc[right, 3] < eff: right -= 1
                y1 = self.group.get_group(each).iloc[right + 1, 3] - eff
                x1 = self.group.get_group(each).iloc[right + 1, 1]
                y2 = self.group.get_group(each).iloc[right, 3] - eff
                x2 = self.group.get_group(each).iloc[right, 1]
                # print(x1,y1,x2,y2)
                mass1 = x1 - y1 * (x2 - x1) / (y2 - y1)
                # print(mass1)
                # print(self.intefuns[each]["fmass"].x)

                left = 0
                while self.group.get_group(each).iloc[left, 3] < eff: left += 1
                y1 = self.group.get_group(each).iloc[left - 1, 3] - eff
                x1 = self.group.get_group(each).iloc[left - 1, 1]
                y2 = self.group.get_group(each).iloc[left, 3] - eff
                x2 = self.group.get_group(each).iloc[left, 1]
                mass2 = x1 - y1 * (x2 - x1) / (y2 - y1)
                # print(mass2)

                # mass1 = newto(lambda x: self.intefuns[each]["feff"](x) - eff,
                #               (max(self.group.get_group(each)[self.map.columns[1]])),h=0.001)
                pi1 = float(self.intefuns[each]["fmass"](mass1))
                # mass2 = newton(lambda x: self.intefuns[each]["feff"](x) - eff,
                #                self.group.get_group(each)[self.map.columns[1]][0])
                pi2 = float(self.intefuns[each]["fmass"](mass2))
                # print(listTodic(self.map.columns, [each, mass1, pi1, eff]))
                newgroup = newgroup.append(listTodic(self.map.columns, [each, mass1, pi1, eff]), ignore_index=True)
                newgroup = newgroup.append(listTodic(self.map.columns, [each, mass2, pi2, eff]), ignore_index=True)

        newgroup = newgroup.sort_values(by=[self.headers[0], self.headers[1]])
        divideg = newgroup.groupby(self.headers[0])
        header = []
        for each, ii in divideg:
            header.append(each)
        g1 = DataFrame(columns=self.map.columns, dtype=float)
        g2 = DataFrame(columns=self.map.columns, dtype=float)
        for jj in header:
            if len(divideg.get_group(jj)) == 2:
                g1 = g1.append(divideg.get_group(jj).iloc[0, :])
                g2 = g2.append(divideg.get_group(jj).iloc[-1, :])
            elif len(divideg.get_group(jj)) < 2:
                g2 = g2.append(divideg.get_group(jj).iloc[-1, :])
        g1.sort_values(by=self.headers[0])
        g2 = g2.sort_values(by=self.headers[0], ascending=False)
        print(g1)
        print(g2)
        self.effmap = self.effmap.append(g1)
        self.effmap = self.effmap.append(g2)
        print(self.effmap)
        # print(newgroup)

    def power(self, speed, pik, Tin=300):
        massflowrate, eff = self.massFlowrate(speed, pik)
        from Engine.GasProperty import cp_Justi, k_Justi
        k = k_Justi(Tin)
        power = massflowrate * cp_Justi(Tin) * Tin * (pow(pik, (k - 1) / k) - 1) / eff
        print("power={} at speed={}, pressure ratio={}, mass flow rate={},efficiency={}".format(power, speed, pik,
                                                                           massflowrate, eff))
        return power

    def pressureRatio(self, speed, power, Tin=300):
        if speed not in self.speeds:
            self.addLine(speed)
        temp = self.group.get_group(speed)[self.headers[2]]

        print("maximum pressure ratio={},minimum pressure ratio={} at speed {}".format(max(temp), min(temp), speed))
        print("\nCalculating the maximum and minimum power of the compressor")
        self.power(speed, min(temp), Tin)
        self.power(speed, max(temp), Tin)
        print("Done!\n")

        PRtemp = (max(temp) + min(temp)) / 2.

        def fun(PR):
            return self.power(speed, PR, Tin) - power

        print("Calculating the pressure ratio")
        h = 0.01
        while abs(fun(PRtemp)) > 1.e-5:
            PRtemp -= fun(PRtemp) / ((fun(PRtemp + h) - fun(PRtemp - h)) / 2 / h)
            print(PRtemp)

        return PRtemp

    # 将插值好的Map写入到文件并打开
    def writedata(self, filename="CompressorMapOut.xlsx", open=True):
        self.map.to_excel(filename)
        if open:
            import os
            os.system("start " + filename)


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
            ax[0] = self.group.get_group(each).plot(x=self.headers[2], y=self.headers[3], label="speed=" + str(each),
                                                    ax=ax[0], marker="x")

        ax[0].set_ylabel(self.headers[3])
        ax[1].set_ylabel(self.headers[1])

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
    

    Turb = Turbine("./TurbineMapMa.xlsx")
    Turb.plot(showlegend=True)
