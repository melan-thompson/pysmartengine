#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/9/16
# @Author  : github.com/Melan_Thompson

from .GasProperty import *
from .Table import ArrayTable,PhsicalVarList
from abc import ABCMeta, ABC
from abc import abstractmethod
from .Algorithm.FileManipulate import get_pakage_dir


def DicplacementVolumn(Bore, Stroke, NumberOfCYlinders=1):
    import math
    return math.pi / 4. * pow(Bore, 2) * Stroke * NumberOfCYlinders


def WibeFunction(theta, SOC, TCD, a=6.908, m=2):
    """Wiebe函数"""
    from math import exp
    if theta < SOC:
        return 0
    else:
        return 1 - exp(-a * pow((theta - SOC) / TCD, m + 1))

def IdealMillerCylcle(T2=400, Rp=0.01, alpha=1.1, eps=18, epsm=15, pik=4, etaTK=0.5, k=None, Hu=42700e3,L0=14.3):
    """
    Thermodynamic Model and Optimization of a Miller Cycle Applied on a Turbocharged Diesel Engine
    T2:进气温度,K
    Rp:预混燃烧比
    alpha:过量空气系数
    eps:压缩比，
    epsm:有效压缩比
    pik:压气机压比，
    etaTK:涡轮增压器总效率
    """
    # from Engine.GasProperty import cv_Justi, cp_Justi, k_Justi
    if k is None:
        k = k_Justi(T2)
    cv = cv_Justi(T2)
    lamda = Rp * Hu / (alpha * L0 * cv * pow(epsm, k - 1) * T2) + 1
    cp = cp_Justi(T2)
    rho = Hu * (1 - Rp) / (k * Hu * Rp + alpha * L0 * pow(epsm, k - 1) * cp * T2) + 1
    pit = pow(1 - pow(eps / epsm, k - 1) * (pow(pik, (k - 1) / k) - 1) / lamda / pow(rho, k) / etaTK, -k / (k - 1))
    eta = 1 - (lamda * pow(rho, k) / pow(eps, k - 1) + (k - 1) * eps / pow(epsm, k) - k / pow(epsm, k - 1) - (k - 1) * (
            1 - pit / pik) * (eps - 1) / pow(epsm, k)) / (lamda - 1 + k * lamda * (rho - 1))
    # BMEP=p0*pik*(((lamda-1)+k*lamda*(rho-1))
    return eta


# def x(theta, start=0, end=60):
#     from math import cos, pi
#     return 0.5 * (1 - cos(pi * (theta - start) / (end - start)))


def FMEP(D: float, cm: float, pme: float, stroke=4) -> float:
    """
    计算平均机械损失压力，适用于四冲程和二冲程
    :param D:Bore，气缸直径(m)
    :param cm:Mean piston moving velocity，活塞平均运动速度(m/s)
    :param pme:Brake mean effective pressure，平均有效压力(Pa)
    :rtype: FMEP，平均机械损失压力(Pa)
    :param stroke:冲程数
    """
    if stroke == 4:
        # return pow(D, -0.1778) * (0.00855 * cm + 0.0789 * pme / 1.e6 - 0.0214) * 1.e6
        return pow(D, -0.1778) * (0.0855 * cm + 0.0789 * pme / 98066.5 - 0.214) * 98066.5
    elif stroke == 2:
        return pow(D, -0.1988) * (0.0829 * cm + 0.08295 * (pme / 98066.5) - 0.3) * 98066.5
    else:
        raise Exception("Stroke number error, either 2 or 4 is allowed")


def FMEP2(cm: float, pme: float) -> float:
    return (0.4 + 0.004 * pme / 1.e5 + 0.09 * cm + 0.0009 * pow(cm, 2)) * 1.e5


def maximumPressure(pme: float):
    """
    由平均有效压力估计最大爆发压力
    :param pme: 平均有效压力，Pa
    :return: 最大爆发压力，Pa
    """
    return (27.65502 + 6.01295 * pme / 1.e5) * 1.e5


def MeEff(D, cm, pme) -> float:
    """
    计算机械效率，适用于四冲程
    :param D:Bore，气缸直径(m)
    :param cm:Mean piston moving velocity，活塞平均运动速度(m/s)
    :param pme:Brake mean effective pressure，平均有效压力(Pa)
    :rtype: Mechanical efficiency,机械效率=BMEP/(BMEP+FMEP)(Pa)
    """
    fmep = FMEP(D, cm, pme)
    print("FMEP={}bar".format(fmep / 1.e5))
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


def BMEP2(IMEP, bore, stroke, speed):
    """
    由IMEP通过机械损失计算BMEP,IMEP*{\eta _m}\left( {{C_m},BMEP} \right) = BMEP
    :param IMEP: 指示平均有效压力,Pa
    :param bore: 缸径,m
    :param stroke: 冲程，m
    :param speed: 转速,rpm
    :return: 平均有效压力
    """
    cm = stroke * speed / 30

    def fun(BMEPpara):
        return BMEPpara - IMEP * MeEff(bore, cm, BMEPpara)

    BMEPresult = IMEP
    h = IMEP / 1.e2
    while abs(fun(BMEPresult) > 1.e-5):
        BMEPresult -= fun(BMEPresult) / ((fun(BMEPresult + h) - fun(BMEPresult - h)) / 2 / h)

    return BMEPresult

def Power(BMEP, n, TVH, stroke=4):
    return BMEP * TVH / (30 * stroke / n)


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
    {m_{f.cycle}} = \frac{{{P_e}{g_f}}}{{3600i}}\frac{{30\tau }}{n}
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


def pimFuel(mfcycle, Vs, phia=1.1, phic=0.85, Tim=300, L0=14.3):
    """
    作用:由燃油喷射量估算进气压力
    公式：{p_{im}} = {m_{f.cycle}}{L_0}\frac{{{\phi _a}}}{{{\phi _c}}}\frac{{{R_g}{T_{im}}}}{{{V_s}}}
    :param mfcycle:Injected mass per cycle per cylinder,单缸循环喷油量(kg)，也可为总的喷油量，此时Vs为发动机总排量
    :param Vs:Single cylinder diceplacement volume,单缸工作容积
    :param phia:Excess air fuel ratio,过量空气系数
    :param phic:Volume metric efficiency充量系数
    :param Tim:进气温度
    :param L0:当量比
    :return:进气压力，Pa
    """
    # from Engine.GasProperty import Rg
    return mfcycle * L0 * phia / phic * Rg(phia) * Tim / Vs


def massAir(Vd, ps=1.e5, Ts=300, VE=1.0, phis=1.0):
    """
    由进气管状态计算发动机进气量
    :param Vd: 排量m^3
    :param ps: 进气管压力(Pa)
    :param Ts: 进气管温度(K)
    :param VE: 充量系数
    :param phis: 扫气系数
    :return: 进气量(kg)
    """
    # from Engine.GasProperty import Rg
    return VE * ps * Vd / Rg() / Ts


def thermalUsageIndex(gi, alpha, Cm, D, S, n, Ps, Ts, strokeNum=4):
    # r'R = \frac{{{{\left( {\alpha {g_i}} \right)}^{0.5}}T_s^{1.5}C_m^{0.78}\left( {0.5 + \frac{D}{{2S}}} \right)}}{{\left( {Dn} \right){{\left( {D{p_s}} \right)}^{0.22}}}}
    # {\xi _T} = \left\{ \begin{array}{l}'
    # 1.028 - 0.00096R,四冲程增压柴油机\\
    # 0.986 - 0.00025R,二冲程直流扫气十字头式增压柴油机
    # \end{array} \right.
    """
    热利用系数
    :param gi:ISFC,指示燃油消耗率(g/(kW*h)),等于机械效率乘以有效燃油消耗率
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
    print("R={}".format(R))
    if strokeNum == 4:
        return 1.028 - 0.00096 * R
    elif strokeNum == 2:
        return 0.986 - 2.5e-4 * R


def simpleWoschini(D, cm, pz, Tz):
    # {K_w} = 0.303{D^{ - 0.214}}{\left( {{v_m}{p_z}} \right)^{0.786}}{T_z}^{ - 0.525}
    """
    简单的Woschi传热公式
    :param D: 缸径(m),
    :param cm:活塞平均速度(m/s)
    :param pz:缸内气体瞬时压力(Pa)
    :param Tz:缸内气体瞬时温度(K)
    :return:传热系数(W/(m^2*K))
    """

    return 0.303 * pow(D, -0.214) * pow(cm * pz / 1.e6, 0.786) * pow(Tz, -0.525)


def chargeTemperature(TBeforeValve=50 + 273.15):
    """
    中速增压四冲程柴油机
    :param TBeforeValve: 进气门前空气的温度(K)
    :return:充量温度(K)
    """
    return 313. + 5. / 6. * (TBeforeValve - 273.15)


def exhaustTemperature(ge=213, etam=0.86, alpha=1.2, phis=1, Ts=300, thermalUsage=0.9, l0=0.495, Hu=42496):
    # {\xi _T}Hu - \frac{{3600}}{{{g_e}{\eta _m}}} + \alpha {\varphi _s}{L_0}{\left( {\mu {c_p}} \right)_s}{T_s} =
    # \left( {\alpha {\varphi _s} - 1 + {\beta _0}} \right){L_0}{\left( {\mu {c_p}} \right)_T}{T_T}
    """
    排气温度预测公式
    :param ge:燃油消耗率
    :param etam:机械效率
    :param alpha:过量空气系数
    :param phis:扫气系数
    :param Ts:进气温度(K)
    :param thermalUsage:
    :param l0:
    :param Hu:
    :return:
    """
    # from Engine.GasProperty import mucps, mucpT
    beta = 1.034
    Tt = 500
    cps = mucps(Ts)

    def fun(Tt):
        return thermalUsage * Hu - 3600 / (ge / 1.e3 * etam) + alpha * phis * l0 * cps * Ts - (
                alpha * phis - 1 + beta) * l0 * mucpT(Tt, alpha, phis) * Tt

    h = 0.1
    while abs(fun(Tt)) > 1.e-5:
        Tt -= fun(Tt) / ((fun(Tt + h) - fun(Tt - h)) / 2 / h)

    return Tt


def exhaustTemperature2(bore, stroke, pme, n, ge, ps, Ts, alpha=1, etav=1, strokeNum=4, l0=0.495, Hu=42496):
    # {\xi _T}Hu - \frac{{3600}}{{{g_e}{\eta _m}}} + \alpha {\varphi _s}{L_0}{\left( {\mu {c_p}} \right)_s}{T_s} = \left( {\alpha {\varphi _s} - 1 + {\beta _0}} \right){L_0}{\left( {\mu {c_p}} \right)_T}{T_T}
    """
    计算排气温度的JDK算法,首先计算机械效率,再计算涡轮前热利用系数，最后迭代求解
    :param bore: 缸径,m
    :param stroke: 冲程,m
    :param pme:平均有效压力,Pa
    :param n: 转速,r/min
    :param ge: 有效燃油消耗率,g/(kW*h)
    :param ps: 进气管压力,Pa
    :param Ts: 进气管温度,K
    :param alpha: 过量空气系数
    :param etav: 扫气系数
    :param strokeNum: 冲程数，默认为4
    :param l0: 燃烧1kg柴油所需的空气量，默认为0.495kmol/kg
    :param Hu: 燃料低热值，kJ/kg
    :return: 排气温度,K
    """
    cm = stroke * n / 30
    print("活塞平均速度{}m/s".format(cm))

    # 计算机械效率
    etam = MeEff(bore, cm, pme)
    print("机械效率{}".format(etam))

    # 计算涡轮前热利用系数
    R = thermalUsageIndex(ge * etam, alpha, cm, bore, stroke, n, ps, Ts, strokeNum)
    print("涡轮前热利用系数{}".format(R))

    # 计算涡前温度
    Tt = exhaustTemperature(ge, etam, alpha, etav, Ts, R, l0, Hu)
    print("涡前温度{}℃".format(Tt - 273.15))

    return Tt


def exhaustTemperature3(bore, stroke, pme, n, ge, Ts, alpha=1, etav=1, effv=0.9, strokeNum=4, l0=0.495, Hu=42496):
    """

    :param bore: 缸径,m
    :param stroke: 冲程,m
    :param pme:平均有效压力,Pa
    :param n: 转速,r/min
    :param ge: 有效燃油消耗率,g/(kW*h)
    :param Ts: 进气管温度,K
    :param alpha: 过量空气系数
    :param etav: 进气管压力,Pa
    :param effv: 扫气系数
    :param strokeNum: 冲程数，默认为4
    :param l0: 燃烧1kg柴油所需的空气量，默认为0.495kmol/kg
    :param Hu: 燃料低热值，kJ/kg
    :return: 排气温度,K
    """

    cm = stroke * n / 30
    print("活塞平均速度{}m/s".format(cm))

    # 计算单缸功率
    VH = 3.1415 / 4 * bore ** 2 * stroke
    power = Power(pme, n, VH, strokeNum)
    print("Power {} W".format(power))

    # 计算喷油量
    mfcylce_temp = mfcycle(ge, power / 1.e3, n)

    ps = pimFuel(mfcylce_temp, VH, alpha, effv, Ts)
    print("进气压力{}bar".format(ps / 1.e5))

    # 计算机械效率
    etam = MeEff(bore, cm, pme)
    print("机械效率{}".format(etam))

    # 计算涡轮前热利用系数
    R = thermalUsageIndex(ge * etam, alpha, cm, bore, stroke, n, ps, Ts, strokeNum)
    print("涡轮前热利用系数{}".format(R))

    # 计算涡前温度
    Tt = exhaustTemperature(ge, etam, alpha, etav, Ts, R, l0, Hu)
    print("涡前温度{}℃".format(Tt - 273.15))

    return Tt


class CylinderGeometry:
    def __init__(self, EngineType=None, bore=None, stroke=None, connecting_rod_length=None, compression_ratio=None,
                 number_of_cylinders=1):
        # 直接由数据库读取
        if EngineType is not None:
            self.EngineType = EngineType
            from pandas import read_excel
            data = read_excel(get_pakage_dir("Engine")+"\\data\\enginedata.xlsx", index_col="机型").loc[EngineType]
            bore = data["bore,mm"] * 1.e-3;
            stroke = data["stroke,mm"] * 1.e-3;
            connecting_rod_length = data["connecting rod length,mm"] * 1.e-3;
            compression_ratio = data["compression ratio"];
            number_of_cylinders = data["number of cylinders"]

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
        self.num_of_cylinders = number_of_cylinders
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
        # 单缸排量
        import math
        return math.pi / 4 * self.bore ** 2 * self.stroke
        # return self.V(180) - self.V(0)

    def totalDisplacedVolume(self):
        return self.num_of_cylinders * self.displacedVolume()

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

        plt.title("Engine:{}\n(bore={}mm,stoke={}mm,rod length={}mm)".format(self.EngineType, self.bore * 1.e3,
                                                                             self.stroke * 1.e3,
                                                                             self.connecting_rod_length * 1.e3),
                  fontsize=10)

        plt.tight_layout()
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


# 放热率抽象类
class HeatRelease(metaclass=ABCMeta):
    def __init__(self):
        self.data = ArrayTable(3, 0)
        self.data.setTableHeader(["Crank angle", "DX", "burned fraction"])
        self.data.setTableUnit(["CA", "1/deg", "/"])

    @abstractmethod
    def DX(self, crank_angle):
        pass

    @abstractmethod
    def X(self, crank_angle):
        pass

    @abstractmethod
    def fractionCA(self, fraction):
        pass

    @abstractmethod
    def plot(self):
        pass

    @abstractmethod
    def updateTable(self, step):
        pass


class SingleWiebe(HeatRelease, ABC):
    def __init__(self, start_of_combstion, total_combustion_duration, m=2., a=6.908):
        super(SingleWiebe, self).__init__()
        assert total_combustion_duration > 0, "Total combustion duration can not less than 0"
        self.SOC = start_of_combstion
        self.TCD = total_combustion_duration
        self.EOC = self.SOC + self.TCD
        self.a = a
        self.m = m

    def updateTable(self, step=0.001):
        self.data.clear()
        from numpy import arange
        for i in arange(self.SOC - 10, self.EOC + 10, step):
            self.data.append([i, self.DX(i), self.X(i)])

    def X(self, crank_angle):
        if crank_angle < self.SOC: return 0
        from math import pow, exp
        temp = -self.a * pow(((crank_angle - self.SOC) / self.TCD), (self.m + 1.))
        return 1 - exp(temp)

    def DX(self, crank_angle):
        if crank_angle < self.SOC: return 0
        import math
        temp = -self.a * pow(((crank_angle - self.SOC) / self.TCD), (self.m + 1.))
        return self.a * (self.m + 1) / self.TCD * pow(((crank_angle - self.SOC) / self.TCD), self.m) * math.exp(
            temp)

    def fractionCA(self, fraction):
        assert 0 < fraction < 1, "burned mass fraction must be between 0 to 1 !"
        return self.data.linearInterpolate(fraction, colx=0, coly=2)

    def plot(self):
        self.updateTable()
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(self.data.table[0].data, self.data.table[1].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()

        plt.subplot(122)
        plt.plot(self.data.table[0].data, self.data.table[2].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()


class DoubleWiebe(HeatRelease, ABC):
    def __init__(self, PF, PSOC, PTCD, DTCD, PM=2, DM=0.8, DSOCin=None):
        super(DoubleWiebe, self).__init__()
        assert 0 <= PF <= 1, "Premixed fraction must be less than 1 and greater than 0, but {} was given".format(PF)
        self.PF = PF
        self.DF = 1 - self.PF
        self.SOC = PSOC
        self.PTCD = PTCD
        self.PM = PM
        self.DTCD = DTCD
        self.DM = DM

        if DSOCin is None:
            self.DSOC = self.SOC + self.PTCD / 2.
        else:
            self.DSOC = DSOCin

        self.EOC = self.DSOC + self.DTCD
        self.PremixWiebe = SingleWiebe(self.SOC, self.PTCD, self.PM)
        self.DisfusionWiebe = SingleWiebe(self.DSOC, self.DTCD, self.DM)

    def updateTable(self, step=0.001):
        self.data.clear()
        from numpy import arange
        for i in arange(self.SOC - 10, self.EOC + 10, step):
            self.data.append([i, self.DX(i), self.X(i)])

    def X(self, CA):
        return self.PF * self.PremixWiebe.X(CA) + self.DF * self.DisfusionWiebe.X(CA)

    def DX(self, CA):
        return self.PF * self.PremixWiebe.DX(CA) + self.DF * self.DisfusionWiebe.DX(CA)

    def fractionCA(self, fraction):
        return self.data.linearInterpolate(fraction, colx=2, coly=0)

    def plot(self):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(self.data.table[0].data, self.data.table[1].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()

        plt.subplot(122)
        plt.plot(self.data.table[0].data, self.data.table[2].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()


# 多Wiebe函数
class MultiWiebe(HeatRelease, ABC):
    def __init__(self, SOC=[-5, 0, 0], TCD=[30, 60, 40], Fraction=[0.1, 0.7, 0.2], m=[2, 3, 5]):
        super(MultiWiebe, self).__init__()
        assert sum(Fraction) == 1, "The sum of fraction must equals to 1"
        self.Fraction = Fraction
        self.SOC = SOC[0]
        self.Wiebe = []
        self.Wiebe.append(SingleWiebe(SOC[0], TCD[0], m[0]))
        SOCtemp = SOC[0]
        for i in range(1, len(SOC)):
            SOCtemp += SOC[i]
            self.Wiebe.append(SingleWiebe(SOCtemp, TCD[i], m[i]))
        self.EOC = max([each.EOC for each in self.Wiebe])
        from numpy import arange
        for i in arange(self.SOC - 10, self.EOC + 10, 0.001):
            self.data.append([i, self.DX(i), self.X(i)])

    def DX(self, crank_angle):
        result = 0.0
        for i in range(len(self.Fraction)):
            result += self.Fraction[i] * self.Wiebe[i].DX(crank_angle)
        return result

    def X(self, crank_angle):
        result = 0.0
        for i in range(len(self.Fraction)):
            result += self.Fraction[i] * self.Wiebe[i].X(crank_angle)
        return result

    def fractionCA(self, fraction):
        return self.data.linearInterpolate(fraction, colx=2, coly=0)

    def plot(self):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(self.data.table[0].data, self.data.table[1].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()
        from numpy import array
        for i in range(len(self.Wiebe)):
            plt.plot(self.Wiebe[i].data.table[0].data, self.Fraction[i] * array(self.Wiebe[i].data.table[1].data),
                     "r-.")

        plt.subplot(122)
        plt.plot(self.data.table[0].data, self.data.table[2].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()


class HeatReleaseData(HeatRelease, ABC):
    def __init__(self, data, SOC=None, EOC=None):
        super(HeatReleaseData, self).__init__()
        # assert data is ArrayTable,"Only ArrayTable class is allowed!!"
        for i in range(data.row):
            self.data.append([data.table[0].data[i], data.table[1].data[i], 0])
        if SOC is None:
            self.SOC = self.__findSOC()
        else:
            self.SOC = SOC

        if EOC is None:
            self.EOC = self.__findEOC()
        else:
            self.EOC = EOC
        print("Start of combustion:{}".format(self.SOC))
        print("End of combustion:{}".format(self.EOC))

        self.data.append([self.SOC, self.data.linearInterpolate(self.SOC, colx=0, coly=1), 0])
        self.data.append([self.EOC, self.data.linearInterpolate(self.EOC, colx=0, coly=1), 0])
        self.data.doQuickSort(0)

        # 把不要的部分切割
        self.data = self.data.slice(self.SOC, self.EOC)
        # self.SOCIndex = self.data.table[0].data.index(self.SOC)
        # self.EOCIndex = self.data.table[0].data.index(self.EOC)

    def updateTable(self, step=0.001):
        pass

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

    def fractionCA(self, fraction):
        assert 0 < fraction < max(self.data.table[2].data), "burn mass fraction out of range!!"
        return self.data.linearInterpolate(fraction, coly=0, colx=2)

    # 归一化
    def nomalization(self):
        # 计算积分
        total = self.data.integrate(1)
        # 将每一个值都做归一化
        for i in range(self.data.row):
            self.data.table[1].data[i] /= total
            self.data.table[2].data[i] = self.data.integrate(1, end_row=i)

    def DX(self, crank_angle):
        return self.data.linearInterpolate(crank_angle, colx=0, coly=1)

    def X(self, crank_angle):
        return self.data.linearInterpolate(crank_angle, colx=0, coly=2)

    # # 积分计算CA10
    # def __analysis(self):
    #     from ArrayTable import PhsicalVarList
    #     burnedfraction = PhsicalVarList([], "Total burned fraction", "/")
    #     self.CA10 = self.CA50 = self.CA90 = None
    #     print("Maximum burned fraction is {}".format(self.data.integrate(1)))
    #     index = 0
    #
    #     while True:
    #         burnedfraction.data.append(self.data.integrate(1, index, 0))
    #         index += 1
    #         if index >= self.data.row:
    #             print("Heat release data has some problem, can not find CA10,CA50,CA90")
    #             return
    #
    #     self.CA10 = self.data.table[0].data[index]
    #     print("Find CA10 at {} °CA".format(self.CA10))
    #
    #     while self.data.integrate(1, index, 0) < 0.5:
    #         index += 1
    #         if index >= self.data.row:
    #             self.data.appendColumn(burnedfraction)
    #             print("Heat release data has some problem, can not find CA50,CA90")
    #             return
    #     self.CA50 = self.data.table[0].data[index]
    #     print("Find CA50 at {} °CA".format(self.CA50))
    #
    #     while self.data.integrate(1, index, 0) < 0.9:
    #         index += 1
    #         if index >= self.data.row:
    #             print("Heat release data has some problem, can not find CA90")
    #             return
    #     self.CA90 = self.data.table[0].data[index]
    #     print("Find CA90 at {} °CA".format(self.CA90))

    def plot(self):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.subplot(121)
        plt.plot(self.data.table[0].data, self.data.table[1].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burn Rate (Normalized by Total Fuel Mass)")
        plt.tight_layout()
        plt.grid()

        plt.subplot(122)
        plt.plot(self.data.table[0].data, self.data.table[2].data)
        plt.xlabel("Crank angle(°CA)")
        plt.ylabel("Burned Fuel (Fraction of Total Fuel Mass)")
        plt.tight_layout()

        plt.grid()
        plt.show()

    def regressWithSingleWiebe(self, plot=True):
        def fun(m):
            result = 0
            wiebe = SingleWiebe(self.SOC, self.EOC - self.SOC, m)
            for i in range(self.data.row):
                result += (wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i]) ** 2
            print("error={}".format(result))
            return result

        # 梯度法
        # it=0
        # m=5;h=0.1
        # while abs(fun(m))>1.e-2 and it<20:
        #     # dfun=((fun(m+h)-fun(m-h))/2./h)
        #     # m-=0.1*dfun
        #     it+=1
        #     m-=fun(m)/((fun(m+h)-fun(m-h))/2./h)

        from sko.PSO import PSO
        from .Algorithm.Optimization import multiopt
        opt=multiopt(fun,lb=[0],ub=[10])
        # opt.GA()
        opt.DE(max_iter=100,size_pop=50)
        # opt.DE()
        x=opt.best_x
        # pso = PSO(func=fun, dim=1, pop=50, max_iter=100, lb=[0], ub=[20])
        # fitness = pso.run()
        print("exp(-{}(x-{})/{})^{}".format(6.908, self.SOC, self.EOC - self.SOC, x[0]))
        print('Absolute error=', x[0])
        print("m={}".format(x[0]))
        result = ArrayTable(3, 0)
        wiebe = SingleWiebe(self.SOC, self.EOC - self.SOC, x[0])
        for i in range(0, self.data.row):
            result.append(
                [self.data.table[0].data[i], self.data.table[1].data[i], wiebe.DX(self.data.table[0].data[i])])
        # result.plot(0, [1, 2])
        if plot:
            result.plot([1, 2])
        return result

    #
    def regressWithDoubleWiebe(self,plot=True):
        def fun(PF, PTCD, DTCD, PM, DM):
            wiebe = DoubleWiebe(PF, self.SOC, PTCD, DTCD, PM, DM, self.EOC - DTCD)
            result = 0
            for i in range(self.data.row):
                result += abs(wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i])
            print("error={}".format(result))
            return result

        # from sko.PSO import PSO
        from .Algorithm.Optimization import multiopt
        lb=[0.0, 10, 10, 0.1, 0.1]
        ub=[1, 90, 60, 5, 10]
        opt=multiopt(fun,lb=lb,ub=ub)
        opt.PSO(50,100)
        opt.GA(50,100)
        opt.PSO(50,100)
        from sko.GA import GA
        # pso = GA(func=fun, n_dim=5, size_pop=100, max_iter=400, lb=[0.0, 10, 10, 0.1, 0.1],
        #           ub=[1, 90, 60, 5, 10])
        # pso = PSO(func=fun, dim=5, pop=100, max_iter=400, lb=[0.0, 10, 10, 0.1, 0.1],
        #           ub=[1, 90, 60, 5, 10], w=0.8, c1=2, c2=2)
        x=opt.best_x

        # fitness = pso.run()
        print("{}*exp(-6.908(x-{})/{})^{}+{}*exp(-6.908(x-{})/{})^{})".format(x[0], self.SOC, x[1],
                                                                              x[3], 1 - x[0],
                                                                              self.EOC - x[2], x[2],
                                                                              x[4]))
        print('best_x is ', x)
        # print('best_y is ', pso.gbest_y)
        result = ArrayTable(3, 0)
        wiebe = DoubleWiebe(x[0], self.SOC, x[1], x[2], x[3], x[4],
                            self.EOC - x[2])
        for i in range(0, self.data.row):
            result.append(
                [self.data.table[0].data[i], self.data.table[1].data[i], wiebe.DX(self.data.table[0].data[i])])
        result.setTableHeader(["Crank angle", "Original HRR", "PSO HRR"])
        result.setTableUnit(["CA", "1/deg", "1/deg"])
        result.plot([1, 2])
        return result
    #
    # def __regressWithSingleWiebeFunc(self, m, a=6.908):
    #     wiebe = SingleWiebe(self.SOC, self.EOC - self.SOC, m, a)
    #     result = 0
    #     for i in range(self.SOCIndex, self.EOCIndex):
    #         result += abs(wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i])
    #     return result
    #
    # def __regressWithDoubleWiebeFunc(self, PF, PTCD, DTCD, PM=2, DM=0.8):
    #     wiebe = DoubleWiebe(PF, self.SOC, PTCD, DTCD, PM, DM, self.EOC - DTCD)
    #     result = 0
    #     for i in range(self.SOCIndex, self.EOCIndex):
    #         result += abs(wiebe.DX(self.data.table[0].data[i]) - self.data.table[1].data[i])
    #     return result

class CylinderPressure:
    def __init__(self, _data, CylGeo=None, move=0):
        self.data = _data
        self.data.setUnitToSI()

        # 转化为以0为上止点曲轴转角
        from .Valve import mod
        for i in range(self.data.row):
            self.data.table[0].data[i] = mod(self.data.table[0].data[i] - move)

        self.data.doQuickSort(0)
        self.CylGeo = CylGeo
        # self.data.plot()

    def FFTTrasform(self, n, tau=4.):
        """
        傅里叶变换分析
        :param n:发动机转速，r/min
        :param tau:发动机冲程
        :return:
        """

        # 计算采样频率
        def basicFrequency(n, CAInterval, stroke=4.):
            t_cycle = 30. * stroke / float(n)
            t_interval = CAInterval * t_cycle / 720.  # 采一次样的时间间隔
            return 1 / t_interval  # 采样频率

        CAInterval = self.data.table[0].data[1] - self.data.table[0].data[0]
        baseF = basicFrequency(n, CAInterval, tau)
        print("Data collection interval is {} CA.".format(CAInterval))
        print("Basic frequency is {} kHz".format(baseF / 1.e3))
        from numpy import fft
        import numpy as np

        dofft = fft.fft(self.data.table[1].data)
        result = ArrayTable(3, 0)
        result.setTableHeader(["Frequency", "Amplitude", "Phase"])
        result.setTableUnit(["Hz", "/", ""])
        result.append([0 * baseF, np.abs(dofft[0]) / self.data.row, np.angle(dofft[0])])
        for i in range(1, self.data.row // 2):
            result.append([i * baseF, np.abs(dofft[i]) / self.data.row * 2, np.angle(dofft[i])])
        return result

    def FFTFilter(self, _cutOffFrequency, amplitudeMutiplier=1):
        from numpy import fft
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

    def CentralMoment(self, start, end, order=0):
        table = self.data.slice(start, end)
        if order == 0:
            return table.integrate()
        if order == 1:
            from .Table import PhsicalVarList
            temp = PhsicalVarList()
            for i in range(table.row):
                temp.append(table.table[0].data[i] * table.table[1].data[i])
            table.appendColumn(temp)
            return table.integrate(2, 0)
        if order > 1:
            thetac = self.CentralMoment(start, end, 1)
            from .Table import PhsicalVarList
            temp = PhsicalVarList()
            for i in range(table.row):
                temp.append(pow(table.table[0].data[i] - thetac, order) * table.table[1].data[i])
            table.appendColumn(temp)
            return table.integrate(2, 0)

    def netHeatReleaseRate(self, CylGeo=None, ivc=None, TBeforeValeve=50 + 173.15, plot=False):
        if CylGeo is None:
            CylGeo = self.CylGeo

        def gamma(T):
            return 1.338 - 6.0e-5 * T + 1.0e-8 * T * T

        if ivc is None:
            ivc = -180 + 50  # 进气门迟后角为30~70

        # from GasProperty import Rg
        mc = self.data.linearInterpolate(ivc, 1) * CylGeo.V(ivc) / Rg() / chargeTemperature(TBeforeValeve)

        result = self.data.selectColumns([0, 1])
        result.diff(1)  # 求dp/dphi
        from .Table import PhsicalVarList
        dQ = PhsicalVarList([], "$dQ/d\varphi$", "J/deg")
        Vdata = PhsicalVarList([], "$V$", "$m^3$")
        Tdata = PhsicalVarList([], "$T$", "K")
        for i in range(result.row):
            p = result.table[1].data[i]

            V = CylGeo.V(result.table[0].data[i]);
            Vdata.data.append(V)

            T = p * V / Rg() / mc;
            Tdata.data.append(T)

            _gamma = gamma(T)

            dQ.data.append(_gamma / (_gamma - 1) * p * CylGeo.DV(result.table[0].data[i]) + 1 / (_gamma - 1) * V *
                           result.table[2].data[i])
        result.appendColumn([Vdata, Tdata, dQ])

        # find start of combustion and end of combustion
        zeroindex = result.findMaxValueIndex(5)
        try:
            while result.table[0].data[zeroindex] > 0:
                zeroindex -= 1
            while result.table[5].data[zeroindex] > 0:
                zeroindex -= 1
            fi0 = result.table[0].data[zeroindex]
            fi1 = result.table[0].data[zeroindex + 1]
            p20 = result.table[5].data[zeroindex]
            p21 = result.table[5].data[zeroindex + 1]
            soc = fi0 - p20 * (fi1 - fi0) / (p21 - p20)
        except:
            print("Can not find start of combustion")
            soc = None

        try:
            zeroindex2 = result.findMaxValueIndex(5)
            while result.table[0].data[zeroindex2] < 0:
                zeroindex2 += 1
            while result.table[5].data[zeroindex2] > 0:
                zeroindex2 += 1
            fi0 = result.table[0].data[zeroindex2 - 1]
            fi1 = result.table[0].data[zeroindex2]
            p20 = result.table[5].data[zeroindex - 1]
            p21 = result.table[5].data[zeroindex2]
            eoc = fi0 - p20 * (fi1 - fi0) / (p21 - p20)
        except:
            print("Can not find end of combustion")
            eoc = None

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(3, 1, figsize=(10, 10))

            ax[0].plot(result.table[0].data, result.table[1].data)
            ax[0].scatter(ivc, result.linearInterpolate(ivc, 1), color="r")
            ax[0].set_ylabel("Cylinder pressure(Pa)")

            ax[1].plot(result.table[0].data, result.table[4].data)

            ax[1].set_ylabel("Temperature(K)")

            ax[2].plot(result.table[0].data, result.table[5].data)
            ax[2].set_ylabel("Net heat release rate(J/deg)")

            plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])

            # ax.axhline(y=0, color='r', linestyle="-.")
            for i in range(3):
                ax[i].axvline(x=0, color='g', linestyle=":")
                ax[i].axvline(x=180, color='g', linestyle=":")
                ax[i].axvline(x=-180, color='g', linestyle=":")
                ax[i].axvline(x=360, color='g', linestyle=":")
                ax[i].axvline(x=-360, color='g', linestyle=":")
            ax[2].axhline(y=0, color='r', linestyle="-.")

            # 最大温度点
            maxTindex = result.findMaxValueIndex(4)
            ax[1].scatter(result.table[0].data[maxTindex], result.table[4].data[maxTindex], color="r")
            ax[1].annotate('Maximum T %.2f K' % result.table[4].data[maxTindex],
                           xy=(result.table[0].data[maxTindex], result.table[4].data[maxTindex]), xycoords='data',
                           xytext=(30, 0), textcoords='offset points',
                           arrowprops=dict(arrowstyle="->"))

            # 燃烧始点和燃烧终点
            if soc is not None and eoc is not None:
                ax[2].scatter(soc, 0, color="r")
                ax[2].annotate('Start of combustion %.3g $^\circ$CA' % soc,
                               xy=(soc, 0), xycoords='data',
                               xytext=(-170, 50), textcoords='offset points',
                               arrowprops=dict(arrowstyle="->"))

                ax[2].scatter(eoc, 0, color="r")
                ax[2].annotate('End of combustion %.4g $^\circ$CA' % eoc,
                               xy=(eoc, 0), xycoords='data',
                               xytext=(0, 50), textcoords='offset points',
                               arrowprops=dict(arrowstyle="->"))

                ax[0].scatter(soc, result.linearInterpolate(soc, 1), color="r")
                ax[0].scatter(eoc, result.linearInterpolate(eoc, 1), color="r")

                ax[1].scatter(soc, result.linearInterpolate(soc, 4), color="r")
                ax[1].scatter(eoc, result.linearInterpolate(eoc, 4), color="r")

            plt.tight_layout()
            plt.show()
            # ax.plot()

        return result, soc, eoc

    def PVDiagram(self, CylGeo=None):
        if CylGeo is None:
            CylGeo = self.CylGeo
        from .Table import ArrayTable
        result = ArrayTable(2, 0)
        for i in range(self.data.row):
            result.append([CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])
        result.setTableHeader(["Cylinder volume", "Cylinder pressure"])
        result.setTableUnit(["$m^3$", "Pa"])
        return result

    # 积分计算IMEP
    def IMEP(self):
        table = self.PVDiagram()
        return table.integrate(1)

    def BMEP(self, speed):
        IMEP = self.IMEP()
        return BMEP2(IMEP,self.CylGeo.bore,self.CylGeo.stroke,speed)

    def MaxPressure(self):
        return max(self.data.table[1].data)

    def Temperature(self, CylGeo):
        from .Table import ArrayTable
        # from GasProperty import Rg
        result = ArrayTable(4, 0)
        for i in range(self.data.row):
            V = CylGeo.V(self.data.table[0].data[i])
            T = self.data.table[1].data[i] * 1.e5 * V / Rg()
            result.append(
                [self.data.table[0].data[i], T, CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])
        result.setTableHeader(["Crank angle", "Cylinder temperature", "Cylinder volume", "Cylinder pressure"])
        result.setTableUnit(["$^\circ CA$", "K", "$m^3$", "bar"])
        return result

    def startOfCombustion(self, type=0, plot=False):
        data = self.data.slice(-60, 60)
        # self.slice(-60, 60)
        data.diff(1)
        data.diff(2)

        # 二阶导数的过零点
        if type == 0:
            zeroindex = data.findMaxValueIndex(3)
            index = data.findMaximumDataIndex(3, 2)
            print(index)
            while data.table[0].data[zeroindex] > 0:
                zeroindex -= 1
            while data.table[3].data[zeroindex] > 0:
                zeroindex -= 1
            fi0 = data.table[0].data[zeroindex]
            fi1 = data.table[0].data[zeroindex + 1]
            p20 = data.table[3].data[zeroindex]
            p21 = data.table[3].data[zeroindex + 1]
            soc = fi0 - p20 * (fi1 - fi0) / (p21 - p20)
            print("Find start of combustion at {} CA".format(soc))

            if plot:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(1, figsize=(10, 10))
                ax.plot(data.table[0].data, data.table[3].data)
                plt.xlabel(data.table[0].ColName + "(" + data.table[0].ColUnit + ")")
                plt.ylabel(data.table[3].ColName + "(" + data.table[3].ColUnit + ")")
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

        # 二阶导数的极大值点
        if type == 1:
            index = self.data.findMaximumDataIndex(3, 1)
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

        # 三阶导数的极大值点
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

    def plot(self, ValveTiming=None):
        if ValveTiming is None:
            from .Valve import ValveDesign
            ValveTiming = ValveDesign()

        self.data.doQuickSort(0)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 10))
        ax.plot(self.data.table[0].data, (self.data.table[1] * 1.e-5).data)

        index = self.data.findMaxValueIndex(1)
        maxpreCA = self.data.table[0].data[index]
        ax.scatter(maxpreCA, self.data.linearInterpolate(maxpreCA, 1))
        ax.annotate('maxium pressure %.5g bar \nat angle %.4g $^\circ$CA' % (
            self.data.linearInterpolate(maxpreCA, 1), maxpreCA),
                    xy=(maxpreCA, self.data.linearInterpolate(maxpreCA, 1)), xycoords='data',
                    xytext=(-0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        plt.xticks([-360, -180, 0, 180, 360], ["-360\nTDC", "-180\nBDC", "0\nTDCF", "180\nBDC", "360\nTDC"])

        # ax.axhline(y=0, color='r', linestyle="-.")
        ax.axvline(x=0, color='g', linestyle=":")
        ax.axvline(x=180, color='g', linestyle=":")
        ax.axvline(x=-180, color='g', linestyle=":")
        ax.axvline(x=360, color='g', linestyle=":")
        ax.axvline(x=-360, color='g', linestyle=":")

        ax.scatter(ValveTiming.IVC, self.data.linearInterpolate(ValveTiming.IVC, 1))
        ax.annotate('IVC %.3g $^\circ$CA' % ValveTiming.IVC,
                    xy=(ValveTiming.IVC, self.data.linearInterpolate(ValveTiming.IVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(ValveTiming.EVO, self.data.linearInterpolate(ValveTiming.EVO, 1))
        ax.annotate('EVO %.3g $^\circ$CA' % ValveTiming.EVO,
                    xy=(ValveTiming.EVO, self.data.linearInterpolate(ValveTiming.EVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(ValveTiming.EVC, self.data.linearInterpolate(ValveTiming.EVC, 1))
        ax.annotate('EVC %.3g $^\circ$CA' % ValveTiming.EVC,
                    xy=(ValveTiming.EVC, self.data.linearInterpolate(ValveTiming.EVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        ax.scatter(ValveTiming.IVO, self.data.linearInterpolate(ValveTiming.IVO, 1))
        ax.annotate('IVO %.3g $^\circ$CA' % ValveTiming.IVO,
                    xy=(ValveTiming.IVO, self.data.linearInterpolate(ValveTiming.IVO, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        soc = self.startOfCombustion()
        ax.scatter(soc, self.data.linearInterpolate(soc, 1))
        ax.annotate('SOC %.3g $^\circ$CA' % soc,
                    xy=(soc, self.data.linearInterpolate(soc, 1)), xycoords='data',
                    xytext=(-100, 0), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

        plt.ylabel("Cylinder pressure(bar)")
        plt.xlabel("Crank angle")
        plt.tight_layout()
        plt.show()

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

    def polyTropicIndex(self, CylGeo, plot=False):
        from .Table import ArrayTable
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

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.plot(result.table[0].data, result.table[1].data)
            ax.axhline(y=-0, color='r', linestyle="-.")
            ax.axhline(y=1, color='r', linestyle="-.")
            ax.axhline(y=1.38, color='r', linestyle="-.")
            # ax.axvline(x=0, color='r', linestyle="-.")
            ax.set_ylim(-5, 5)
            ax.set_xlim(-120, 160)

            plt.yticks([0, 1, 1.4], ["n=0,dp=0", "n=1,dT=0", "n=$\gamma$,dQ=0"])

            plt.tight_layout()
            plt.show()

        return result

    def LogP_LogV(self, CylGeo=None, ivc=-180 + 50, evo=180 - 10, phline=None, plot=False):
        if CylGeo is None:
            CylGeo = self.CylGeo

        heatRelease, soc, eoc = self.netHeatReleaseRate(CylGeo, ivc)
        # soc=self.startOfCombustion(type=0)

        from math import log
        result = self.data.createSimilarEmptyTable([0, 1])
        for i in range(self.data.row):
            result.append([log(CylGeo.V(self.data.table[0].data[i])), log(self.data.table[1].data[i])])

        compressPFi = self.data.slice(ivc, soc)
        compress = self.data.createSimilarEmptyTable([0, 1])
        for i in range(compressPFi.row):
            compress.append([log(CylGeo.V(compressPFi.table[0].data[i])), log(compressPFi.table[1].data[i])])

        # 拟合
        compressfit, indexcom = compress.polyFit(1)

        intakeProcess = self.data.slice(-360, -180)
        # intakeProcess.show()
        # intakeProcess.plot()
        # intakeProcess.show()
        intake = self.data.createSimilarEmptyTable([0, 1])
        for i in range(intakeProcess.row * 9 // 13):
            intake.append([log(CylGeo.V(intakeProcess.table[0].data[i])), log(intakeProcess.table[1].data[i])])

        # 拟合
        intakefit, indexintake = intake.polyFit(1, _deg=0)

        expensePFi = self.data.slice(eoc, evo)
        expense = self.data.createSimilarEmptyTable([0, 1])
        for i in range(expensePFi.row):
            expense.append([log(CylGeo.V(expensePFi.table[0].data[i])), log(expensePFi.table[1].data[i])])

        expensefit, indexex = expense.polyFit(1)

        # compressfit.plot(2)

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, figsize=(20, 10))
            ax.plot(result.table[0].data, result.table[1].data)

            ax.scatter(log(CylGeo.V(soc)), log(self.data.linearInterpolate(soc, 1)), color="g")
            ax.annotate('SOC %.2f $^\circ$CA' % soc,
                        xy=(log(CylGeo.V(soc)), log(self.data.linearInterpolate(soc, 1))),
                        xycoords='data',
                        xytext=(-40, -50), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.scatter(log(CylGeo.V(eoc)), log(self.data.linearInterpolate(eoc, 1)), color="g")
            ax.annotate('EOC %.2f $^\circ$CA' % eoc,
                        xy=(log(CylGeo.V(eoc)), log(self.data.linearInterpolate(eoc, 1))),
                        xycoords='data',
                        xytext=(0, 20), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.scatter(log(CylGeo.V(ivc)), log(self.data.linearInterpolate(ivc, 1)), color="g")
            ax.annotate('IVC %.1f $^\circ$CA' % ivc,
                        xy=(log(CylGeo.V(ivc)), log(self.data.linearInterpolate(ivc, 1))),
                        xycoords='data',
                        xytext=(-40, 40), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.scatter(log(CylGeo.V(evo)), log(self.data.linearInterpolate(evo, 1)), color="g")
            ax.annotate('EVO %.1f $^\circ$CA' % evo,
                        xy=(log(CylGeo.V(evo)), log(self.data.linearInterpolate(evo, 1))),
                        xycoords='data',
                        xytext=(-100, 20), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            # 压缩过程指数
            ax.plot(compressfit.table[0].data, compressfit.table[2].data, "r-.")
            # compressfit.plot()
            compressfit.doQuickSort()
            ax.annotate('compress isentropic index %.2f' % -indexcom[0],
                        xy=((compressfit.table[0].data[0] + compressfit.table[0].data[-1]) / 2,
                            compressfit.linearInterpolate(
                                (compressfit.table[0].data[0] + compressfit.table[0].data[-1]) / 2, 2)),
                        xycoords='data',
                        xytext=(-150, -50), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            # 膨胀过程指数
            ax.plot(expensefit.table[0].data, expensefit.table[2].data, "r-.")
            ax.annotate('expense isentropic index %.2f' % -indexex[0],
                        xy=((expensefit.table[0].data[0] + expensefit.table[0].data[-1]) / 2,
                            expensefit.linearInterpolate(
                                (expensefit.table[0].data[0] + expensefit.table[0].data[-1]) / 2, 2)),
                        xycoords='data',
                        xytext=(20, 10), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            ax.plot(intakefit.table[0].data, intakefit.table[2].data, "r-.")
            from math import exp
            ax.annotate('pressure %.2f bar' % (exp(indexintake[0]) / 1.e5),
                        xy=((intakefit.table[0].data[0] + intakefit.table[0].data[-1]) / 2,
                            intakefit.linearInterpolate(
                                (intakefit.table[0].data[0] + intakefit.table[0].data[-1]) / 2, 2)),
                        xycoords='data',
                        xytext=(-20, -20), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))

            plt.xlabel("$log V$")
            plt.ylabel("$log p$")

            # if phline is not None:
            #     ax.axhline(y=log(phline), color='r', linestyle="-.")

            plt.tight_layout()
            plt.show()

        return result

    def EquivalentIsotropicIndex(self, CylGeo, reference_point=-200):
        from .Table import ArrayTable
        from math import log
        result = ArrayTable(2, 0)
        Vref = CylGeo.V(reference_point)
        pref = self.data.linearInterpolate(reference_point, 1)
        result.setTableHeader(["Crank angle", "Equivalent Isotropic Index"])
        result.setTableUnit(["°CA", "/"])
        for i in range(self.data.row - 1):
            temp = log(self.data.table[1].data[i + 1] / pref) / log(Vref / CylGeo.V(self.data.table[0].data[i + 1]))
            result.append([self.data.table[0].data[i], temp])

        return result

    def heatReleaseRate(self):
        pass

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
            from .Table import PhsicalVarList
            aftersmoothcol = PhsicalVarList()
            aftersmoothcol.ColName = "Pressure after smooth"
            aftersmoothcol.ColUnit = self.data.table[1].ColUnit
            aftersmoothcol.data = aftersmooth
            self.data.table.append(aftersmoothcol)
            self.data.col += 1
            return self.data



class MillerCycle:
    def __init__(self, CylGeo, p0=101325, T0=273.15 + 20):
        self.p0 = p0
        self.T0 = T0
        self.CylGeo = CylGeo
        
        table = ArrayTable(5, 0)
        table.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        table.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])
        self.data = table

    def initCompressor(self, etak, pik):
        self.pik = pik
        self.etak = etak
        self.pk = self.p0 * pik
        from Engine.Compressor import TAfterCompressor
        self.Tk = TAfterCompressor(self.T0, pik, etak)

    def intake(self, IVC=-30, phic=1):
        # 这里IVC为进气门相对于下止点的关闭角度，为0时代表在下止点关闭
        self.IVC = IVC + 180
        print("Effective compression ratio is {}".format(self.CylGeo.V(180 + IVC) / self.CylGeo.V(0)))
        if IVC <= -90 or IVC > 150:
            raise Exception("Intake valve close timing can not earlier than -90 and later than 150 with respect to "
                            "BDC,but {} is given".format(IVC))
        from Engine.GasProperty import Rg
        from numpy import arange
        for i in arange(0, 180 + IVC, 1):
            V = self.CylGeo.V(i)
            self.data.append([i, self.CylGeo.V(i), self.pk, self.Tk, phic * self.pk * V / Rg(1.e8) / self.Tk])

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
        from Engine.GasProperty import DieselMixture
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
        from Engine.GasProperty import DieselMixture
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
        from Engine.Compressor import piT
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
            # print(piTinit)

        # h = 0.01
        # while abs(fun(piTinit)) > 1.e-2:
        #     dfun = (fun(piTinit + h) - fun(piTinit - h)) / 2 / h
        #     piTinit -= fun(piTinit) / dfun
        #     print(piTinit)
        self.pit = piTinit
        self.T7 = self.p0 * piTinit / self.p6 * self.T6
        self.p7 = self.p0 * piTinit
        self.data.append(
            [540., self.CylGeo.V(540), self.p7, self.T7, self.p7 * self.CylGeo.V(540) / Rg(self.alpha) / self.T7])

    def subcritical(self):
        from numpy import arange
        for i in arange(540 + 0.1, 720, 1):
            V = self.CylGeo.V(i)
            self.data.append([i, V, self.p7, self.T7, self.p7 * V / Rg(self.alpha) / self.T7])

    def analyze(self, plot=True, speed=2000, tau=4, Hu=42700e3, index=0):
        work = self.data.integrate(_coly=2, _colx=1)
        effi = work / self.mix.gf / Hu
        power = work * self.CylGeo.num_of_cylinders / (30 * tau / speed)
        IMEP = work / self.CylGeo.displacedVolume()
        print("Indicated power={}kW".format(power / 1.e3))
        print("IMEP={} bar".format(IMEP / 1.e5))
        print("Indicated thermal efficiency={}".format(effi))

        if plot:
            self.plot(index)
        return effi, IMEP, power

    def plot(self, index=1):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 9))

        temp = [each / 1.e5 for each in self.data.table[2].data]

        plt.rcParams['font.sans-serif'] = ['SimHei']

        ax[0].plot(self.data.table[index].data, temp)
        ax[0].axhline(y=self.p0 / 1.e5, color='r', linestyle="-.")
        ax[0].set_xlabel("曲轴转角℃A")
        ax[0].set_ylabel("气缸压力,bar")
        ax[1].set_xlabel("曲轴转角℃A")
        ax[1].set_ylabel("缸内温度,℃")
        ax[2].set_xlabel("曲轴转角℃A")
        ax[2].set_ylabel("缸内质量,kg")

        ax[1].plot(self.data.table[index].data, self.data.table[3].data)
        ax[1].axhline(y=self.T0, color='r', linestyle="-.")

        ax[2].plot(self.data.table[index].data, self.data.table[4].data)
        # ax[2].axhline(y=self.T0, color='r', linestyle="-.")

        style = dict(size=10, color='black')
        ax[0].text(self.data.table[1].data[0], -self.p0 / 1.e5, "$p_0=%.2fbar$" % (self.p0 / 1.e5), ha="left", **style)

        plt.tight_layout()
        plt.show()



class HeatTransfer:
    def __init__(self, Head_Temperature, Piston_Temperature, Cylinder_Temperature, CorrectCoeff=1):
        self.HT = Head_Temperature
        self.PT = Piston_Temperature
        self.CT = Cylinder_Temperature
        self.Coeff = CorrectCoeff

    def heatTransferRate(self, D, S, V, N, p, T):
        from math import pi, pow
        cm = N * S / 30
        AL = 4. * V / D
        # print(simpleWoschini(D, cm, p, T))
        return self.Coeff * simpleWoschini(D, cm, p, T) * (
                AL * (T - self.CT) + pi * pow(D, 2) / 4. * (2 * T - self.HT - self.PT)) / (6 * N)


from .GasProperty import Rg, k_Justi, cp_Justi, cv_Justi, u_Justi, h_Justi, m_fuel



class Volume:
    def __init__(self, V=1, p0=1.e5, T0=20+273.15, alphak=1.e8):
        # initialize the volume state,初始化容积、质量、温度、压力
        self.V = V
        self.m = p0 * V / T0 / Rg(alphak)
        self.T = T0
        self.p = p0
        self.alphak = alphak

        # 初始化Rg,k,cp,cv,u,h,a
        self.Rg = Rg(self.alphak)
        self.k = k_Justi(self.T, self.alphak)
        self.cp = cp_Justi(self.T, self.alphak)
        self.cv = cv_Justi(self.T, self.alphak)
        self.u = u_Justi(self.T, self.alphak)
        self.h = h_Justi(self.T, self.alphak)
        self.a=pow(self.k*self.Rg*self.T,0.5)

        # 记录下这一步的质量变化
        self.dm_record = None

        # 连接
        self.last = None
        self.next = None

        self.data = ArrayTable(11, 0)
        self.data.setTableHeader(["时间", "质量", "压力", "温度", "过量空气系数", "定压比热", "定容比热", "气体常数", "比热容比", "热力学能", "焓"])

    def dm(self):
        if self.last is None:
            min = 0
        else:
            min = self.last.dm_record
        if self.next is None:
            mout = 0
        else:
            mout = self.next.dm_record
        self.dm_record = min - mout
        return self.dm_record

    # 必须先算dm再算dp
    def dp(self, n=1.4):
        return self.k * Rg(self.alphak) * self.T / self.V * self.dm_record

    def dT(self):
        if self.last is None:
            Hin = 0
        elif self.last.dm_record > 0:
            Hin = self.last.last.h * self.last.dm_record
        elif self.last.dm_record <= 0:
            Hin = self.h * self.last.dm_record

        if self.next is None:
            Hout = 0
        elif self.next.dm_record > 0:
            Hout = self.h * self.next.dm_record
        elif self.next.dm_record <= 0:
            Hout = self.next.next.h * self.next.dm_record

        return 1 / self.m / self.cv * (Hin - Hout - self.u * self.dm_record)

    def dAlpha(self, L0=14.3):
        if self.last is None:
            alphaIn = 0
            dmin = 0
        elif self.last.dm_record > 0:
            alphaIn = self.last.last.alphak
            dmin = self.last.dm_record
        elif self.last.dm_record <= 0:
            alphaIn = self.alphak
            dmin = self.last.dm_record

        if self.next is None:
            alphaOut = 0
            dmout = 0
        elif self.next.dm_record > 0:
            alphaOut = self.alphak
            dmout = self.next.dm_record
        elif self.next.dm_record <= 0:
            alphaOut = self.next.next.alphak
            dmout = self.next.dm_record

        temp1 = (self.alphak * L0 + 1) / (alphaIn * L0 + 1)
        temp2 = (self.alphak * L0 + 1) / (alphaOut * L0 + 1)
        return (dmin * (1 - temp1) - dmout * (1 - temp2)) / (m_fuel(self.m, self.alphak, L0) * L0)

    def update(self):
        # 温度压力变化后容积的变化
        self.Rg = Rg(self.alphak)
        self.k = k_Justi(self.T, self.alphak)
        self.cp = cp_Justi(self.T, self.alphak)
        self.cv = cv_Justi(self.T, self.alphak)
        self.u = u_Justi(self.T, self.alphak)
        self.h = h_Justi(self.T, self.alphak)
        self.p = self.m * self.Rg * self.T / self.V
        self.a=pow(self.k*self.Rg*self.T,0.5)

    def record(self, t):
        self.update()
        self.data.append([t, self.m, self.p, self.T, self.alphak, self.cp, self.cv, self.Rg, self.k, self.u, self.h])


class Inj:
    def __init__(self, mf, L0=14.3, QH=42700e3):
        self.mf = mf
        self.L0 = L0
        self.QH = QH


class EngCylinder(Volume):
    def __init__(self, Fi, pin, Tin, AFAin, speed, cylgeo, hrr, cylht, inj):
        self.Fi = Fi
        self.p = pin
        self.T = Tin
        self.alphak = AFAin
        self.N = speed

        self.Geo = cylgeo
        self.HRR = hrr
        self.CylHT = cylht
        self.inj = inj

        self.m = self.p * self.Geo.V(Fi) / Rg(AFAin) / self.T
        # self.Xk = self.m / (1 + self.inj.L0 * self.alphak) / self.inj.mf

        self.update()
        # self.DX = self.HRR.DX(self.Fi)

        # 记录下这一步的质量变化
        self.dm_record = None

        # 连接
        self.last = None
        self.next = None

        self.data = ArrayTable(19, 0)
        self.data.setTableHeader(["曲轴转角", "气缸容积", "缸内压力", "缸内温度", "广义过量空气系数",
                                  "缸内质量", "废气质量", "空气质量", "进入质量流量", "出质量流量",
                                  "定容比热", "定压比热", "气体常数", "比热容比",
                                  "燃烧放热率", "气缸传热率", "进气焓流", "排气焓流", "对外做功"])
        self.data.setTableUnit(["CA", "m^3", "Pa", "K", "/",
                                "kg", "kg", "kg", "kg/s", "kg/s",
                                "J/(kg*K)", "J/(kg*K)", "J/(kg*K)", "/",
                                "W", "W", "W", "W", "W"])

    # 计算缸内质量变化
    def dm(self):
        if self.last is None:
            self.__dmin = 0
        else:
            self.__dmin = self.last.dm_record

        if self.next is None:
            self.__dmout = 0
        else:
            self.__dmout = self.next.dm_record
        self.__DX = self.HRR.data.linearInterpolate(self.Fi, 1, 1)
        self.dm_record = self.__dmin - self.__dmout + self.inj.mf * self.__DX
        return self.dm_record

    # 计算缸内成分变化
    def dAlpha(self):
        # self.Xk = self.m / (1 + self.inj.L0 * self.alphak) / self.inj.mf

        if self.last is None:
            alphaIn = self.alphak
        elif self.last.dm_record > 0:
            alphaIn = self.last.last.alphak
        elif self.last.dm_record <= 0:
            alphaIn = self.alphak

        if self.next is None:
            alphaOut = self.alphak
        elif self.next.dm_record > 0:
            alphaOut = self.alphak
        elif self.next.dm_record <= 0:
            alphaOut = self.next.next.alphak

        temp1 = (self.alphak * self.inj.L0 + 1) / (alphaIn * self.inj.L0 + 1)
        temp2 = (self.alphak * self.inj.L0 + 1) / (alphaOut * self.inj.L0 + 1)
        return (self.__dmin * (1 - temp1) - self.__dmout * (1 - temp2) - self.__DX * (
                self.inj.L0 * self.alphak + 1) * self.inj.mf) / (m_fuel(self.m, self.alphak, self.inj.L0) * self.inj.L0)

    def dT(self):
        # self.update()
        if self.last is None:
            self.__Hin = 0
        elif self.last.dm_record > 0:
            self.__Hin = self.last.last.h * self.last.dm_record
        elif self.last.dm_record <= 0:
            self.__Hin = self.h * self.last.dm_record

        if self.next is None:
            self.__Hout = 0
        elif self.next.dm_record > 0:
            self.__Hout = self.h * self.next.dm_record
        elif self.next.dm_record <= 0:
            self.__Hout = self.next.next.h * self.next.dm_record

        self.__DW = self.p * self.Geo.DV(self.Fi)
        self.__DQ = self.__DX * self.inj.mf * self.inj.QH
        self.__DQW = self.CylHT.heatTransferRate(self.Geo.bore, self.Geo.stroke, self.V, self.N, self.p, self.T)

        return (
                       self.__DQ + self.__Hin - self.__Hout - self.__DQW - self.__DW - self.u * self.dm_record) / self.cv / self.m

    def update(self):
        # 温度压力变化后容积的变化
        self.V = self.Geo.V(self.Fi)
        self.Rg = Rg(self.alphak)
        self.k = k_Justi(self.T, self.alphak)
        self.cp = cp_Justi(self.T, self.alphak)
        self.cv = cv_Justi(self.T, self.alphak)
        self.u = u_Justi(self.T, self.alphak)
        self.h = h_Justi(self.T, self.alphak)
        self.p = self.m * self.Rg * self.T / self.V

    def M_exhaust(self):
        return self.m * (1 + self.inj.L0) / (self.inj.L0 * self.alphak + 1)

    def M_air(self):
        return self.m - self.M_exhaust()

    def RecordThisStep(self):
        self.update()
        self.data.append(
            [self.Fi, self.V, self.p, self.T, self.alphak, self.m, self.M_exhaust(), self.M_air(), self.__dmin,
             self.__dmout, self.cv, self.cp, self.Rg, self.k, self.__DQ, self.__DQW, self.__Hin, self.__Hout,
             self.__DW])
        return True

    def plot(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(5, 3, figsize=(20, 10))
        ax[0][0].plot(self.data.table[0].data, self.data.table[2].data)
        ax[0][0].set_ylabel(self.data.table[2].ColName + "," + self.data.table[2].ColUnit)

        ax[0][1].plot(self.data.table[0].data, self.data.table[3].data)
        ax[0][1].set_ylabel(self.data.table[3].ColName + "," + self.data.table[3].ColUnit)

        ax[0][2].plot(self.data.table[0].data, self.data.table[4].data)
        ax[0][2].set_ylabel(self.data.table[4].ColName + "," + self.data.table[4].ColUnit)

        # 质量
        ax[1][0].plot(self.data.table[0].data, self.data.table[5].data)
        ax[1][0].set_ylabel(self.data.table[5].ColName + "," + self.data.table[5].ColUnit)

        ax[1][1].plot(self.data.table[0].data, self.data.table[6].data)
        ax[1][1].set_ylabel(self.data.table[6].ColName + "," + self.data.table[6].ColUnit)

        ax[1][2].plot(self.data.table[0].data, self.data.table[7].data)
        ax[1][2].set_ylabel(self.data.table[7].ColName + "," + self.data.table[7].ColUnit)

        # 能量
        ax[2][0].plot(self.data.table[0].data, self.data.table[14].data)
        ax[2][0].set_ylabel(self.data.table[14].ColName + "," + self.data.table[14].ColUnit)

        ax[2][1].plot(self.data.table[0].data, self.data.table[15].data)
        ax[2][1].set_ylabel(self.data.table[15].ColName + "," + self.data.table[15].ColUnit)

        ax[2][2].plot(self.data.table[0].data, self.data.table[16].data)
        ax[2][2].set_ylabel(self.data.table[16].ColName + "," + self.data.table[16].ColUnit)

        ax[3][0].plot(self.data.table[0].data, self.data.table[17].data)
        ax[3][0].set_ylabel(self.data.table[17].ColName + "," + self.data.table[17].ColUnit)

        ax[3][1].plot(self.data.table[0].data, self.data.table[18].data)
        ax[3][1].set_ylabel(self.data.table[18].ColName + "," + self.data.table[18].ColUnit)

        ax[3][2].plot(self.data.table[0].data, self.data.table[9].data)
        ax[3][2].set_ylabel(self.data.table[9].ColName + "," + self.data.table[9].ColUnit)

        ax[4][0].plot(self.data.table[0].data, self.data.table[8].data)
        ax[4][0].set_ylabel(self.data.table[8].ColName + "," + self.data.table[8].ColUnit)

        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 显式负号
        plt.rcParams['axes.unicode_minus'] = False
        plt.tight_layout()
        plt.show()
        # self.data.plot(2)


# 计算能量平衡
class EnergyBalance:
    def __init__(self, bore, stroke, n, ge, Pe, Ts, ps, alpha=2, NOC=1, etav=1.4,piT=2,etaT=0.6,T0=20+273.15,Hu=42496, L0=14.3, linkType=1):
        """
        :param bore:缸径,m
        :param stroke:冲程,m
        :param n: 发动机转速,r/min
        :param ge:发动机油耗,g/(kW*h)
        :param Pe:发动机功率，kW
        :param Ts:进气温度,K
        :param ps:进气压力,Pa
        :param alpha:过量空气系数
        :param NOC:气缸数
        :param etav:扫气系数
        :param Hu:燃料低热值,kJ/kg
        :param L0:燃空当量比,柴油14.3，汽油14.7
        """
        self.nodes = [{"name": "燃料功率"},
                      {"name": "指示功"}, {"name": "曲轴输出功"}, {"name": "机械损失"}, {"name": "气缸排气能量"}, {"name": "气缸冷却水损失"}
            , {"name": "余项损失"},{"name": "废气损失"},{"name": "中冷损失"}]

        self.Pe = Pe
        self.VH = 3.1415 / 4 * bore ** 2 * stroke * NOC
        print("排量{}L".format(self.VH * 1.e3))

        # 计算平均有效压力
        self.pme = Pe * 1.e3 * 120 / n / self.VH
        print("平均有效压力:{}bar".format(self.pme / 1.e5))

        # 计算活塞平均速度
        self.cm = stroke * n / 30
        print("活塞平均速度：{}".format(self.cm))

        # 计算机械效率
        self.etam = MeEff(bore, self.cm, self.pme)
        print("机械效率:{}".format(self.etam))

        # 计算指示功率
        self.Pi = self.Pe / self.etam

        # 计算机械损失功率
        self.Pm = self.Pi - self.Pe

        self.mf = ge * Pe / 3600 / 1.e3  # kg/s
        self.ma = self.mf * alpha * L0  # kg/s
        self.me = self.mf + self.ma  # kg/s

        self.Ha = cp_Justi(Ts) * self.ma * Ts / 1.e3

        self.Pf = self.mf * Hu

        # 计算涡轮前热利用系数
        R = thermalUsageIndex(ge * self.etam, alpha, self.cm, bore, stroke, n, ps, Ts, 4)
        print("涡轮前热利用系数{}".format(R))

        # 计算涡前温度
        Tt = exhaustTemperature(ge, self.etam, alpha, etav, Ts, R)
        print("涡前温度{}℃".format(Tt - 273.15))
        He = cp_Justi(Tt, alpha) * self.me * Tt / 1.e3

        Table = ArrayTable()
        Table.readCSVFile("EnergyBalance.csv")
        qw = Table.GPR([bore * 1.e3, n, self.pme / 1.e5], [1, 4, 5], 9, trained=False)[0][0] / 1.e2
        print(qw)

        from .Compressor import TAfterTurbine
        Tt0=TAfterTurbine(Tt,piT,etaT)
        print("涡轮后温度{}℃".format(Tt0-273.15))
        He0=cp_Justi(Tt, alpha) * self.me * Tt0 / 1.e3

        H0=cp_Justi(T0) * self.ma * Ts / 1.e3

        Qincool=He-self.Ha-(He0-H0)

        if linkType == 0:
            self.links = [
                # {"source": self.nodes[0]["name"], "target": self.nodes[2]["name"], "value": round(self.Pf, 4)},
                # {"source": self.nodes[0]["name"], "target": self.nodes[2]["name"], "value": round(self.Ha, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[3]["name"], "value": round(self.Pi, 4)},
                {"source": self.nodes[3]["name"], "target": self.nodes[6]["name"], "value": round(self.Pm, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[4]["name"], "value": round(He - self.Ha, 4)},
                {"source": self.nodes[3]["name"], "target": self.nodes[7]["name"], "value": round(self.Pe, 4)}

            ]
        else:
            self.links = [
                # {"source": self.nodes[0]["name"], "target": self.nodes[2]["name"], "value": round(self.Pf, 4)},
                # {"source": self.nodes[0]["name"], "target": self.nodes[2]["name"], "value": round(self.Ha, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[1]["name"],
                 "value": round(self.Pi / self.Pf, 4)},
                {"source": self.nodes[1]["name"], "target": self.nodes[3]["name"],
                 "value": round(self.Pm / self.Pf, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[4]["name"],
                 "value": round((He - self.Ha) / self.Pf, 4)},
                {"source": self.nodes[1]["name"], "target": self.nodes[2]["name"],
                 "value": round(self.Pe / self.Pf, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[5]["name"], "value": round(qw, 4)},
                {"source": self.nodes[0]["name"], "target": self.nodes[6]["name"],
                 "value": round(1 - self.Pi / self.Pf - (He - self.Ha) / self.Pf - round(qw, 4), 4)},
                {"source": self.nodes[4]["name"], "target": self.nodes[8]["name"], "value": round(Qincool/self.Pf, 4)},
                {"source": self.nodes[4]["name"], "target": self.nodes[7]["name"], "value": round((He0-H0)/self.Pf, 4)},

            ]

    def plot(self, open=True):
        from pyecharts.charts import Sankey
        from pyecharts import options as opts

        pic = (
            Sankey().add('', self.nodes, self.links,
                         linestyle_opt=opts.LineStyleOpts(opacity=0.3, curve=0.5, color="source"),
                         label_opts=opts.LabelOpts(position="top"), node_gap=30, ).set_global_opts(
                title_opts=opts.TitleOpts(title='')))

        pic.render('test2.html')
        if open == True:
            import os
            os.system("start test2.html")




# 插值函数
def xi(theta, theta0, theta1):
    if theta < theta0 or theta > theta1:
        raise Exception("Wrong error!!")
    from math import cos, pi
    return 0.5 * (1 - cos(pi * (theta - theta0) / (theta1 - theta0)))


class CylPreReconstruct:
    """
    缸压重构的过程，
    self.CylGeo:气缸几何参数
    self.ValveTiming:进排气门正时

    self.mix:缸内的混合气体

    self.T2:燃烧始点的温度
    self.p2:燃烧始点的压力
    """

    def __init__(self, CylGeo, ValveTiming=None, p0=1.e5, T0=20 + 273.15, pe0=None):
        """
        初始化气缸几何参数和气阀参数
        :param CylGeo: 气缸几何参数
        :param ValveTiming: 进排气正时参数，如果没有则为默认值
        :param p0:环境压力，Pa
        :param T0:环境温度，K
        :param pe0:排气压力，正常情况下等于进气压力，考虑到潜艇等有被压存在
        """
        from .Valve import ValveDesign
        self.CylGeo = CylGeo

        # 气门正时
        if ValveTiming is None:
            self.ValveTiming = ValveDesign()
        else:
            self.ValveTiming = ValveTiming

        self.p0 = p0
        self.T0 = T0

        # 排气被压
        if pe0 is None:
            self.pe0 = p0
        else:
            self.pe0 = pe0

        # 记录压缩过程
        self.CompressData = ArrayTable(5, 0)
        self.CompressData.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        self.CompressData.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])

        # 记录膨胀过程
        self.ExpanseData = ArrayTable(5, 0)
        self.ExpanseData.setTableHeader(["Crank angle", "V", "p", "T", "m"])
        self.ExpanseData.setTableUnit(["°CA", "m^3", "Pa", "K", "kg"])

        # 重构的缸压
        self.Rpressure = ArrayTable(2, 0)

        # 换气过程
        self.GasExchange = ArrayTable(2, 0)

        # 最终的缸压
        self.data = ArrayTable(2, 0)
        self.data.setTableHeader(["Crank angle", "Cylinder pressure"])
        self.data.setTableUnit(["°CA", "Pa"])

        # 排气温度
        self.Tt = None

    # 由充量系数和过量空气系数计算燃油消耗率，首先假设一个Ps
    def compress(self, pim=2.5e5, Tim=300, Tr=300, xr=0.0, kc=1.3, phic=1):
        """
        计算压缩过程，即从进气门关闭到燃烧开始时缸内的过程，需要估算排气温度
        :param pim:进气管压力，Pa
        :param Tim: 进气管温度,K
        :param Tr: 排气温度，K
        :param xr: 废气系数
        :param kc: 压缩多变指数
        :param phic: 充量系数
        :return:
        """
        from .GasProperty import DieselMixture, Rg
        self.kc = kc
        ivc = self.ValveTiming.IVC
        self.pim = pivc = pim

        Tim = 313. + 5. / 6. * (Tim - 273.15)  # 新鲜充量温度估算经验公式
        self.Tim = Tivc = Tim * (1 - xr) + xr * Tr

        # 单缸的新鲜充量质量
        mivc = phic * pivc * self.CylGeo.V(180) / Rg() / Tivc
        self.mix = DieselMixture()

        # 由进气量、残余废气系数和循环喷油量初始化缸内气体，先随意给定
        self.mix.init_With_Mc_r(mivc, xr, (1 - xr) * mivc / 14.3 / 2)
        print("Intake air mass {} mg".format(mivc * 1.e6))
        print("Tivc=", Tivc)

        Vivc = self.CylGeo.V(ivc)

        # 压缩过程线
        from Valve import mod
        from numpy import arange
        for i in arange(ivc, ivc + 720):
            V = self.CylGeo.V(i)
            T = Tivc * pow(Vivc / V, kc - 1)
            p = pivc * pow(Vivc / V, kc)
            m = p * V / self.mix.Rg_gas(0) / T
            self.CompressData.append([mod(i), V, p, T, m])
        self.CompressData.doQuickSort(0)

        # for i in range(self.CompressData.row):
        #     self.CompressData.table[0].data[i] = mod(self.CompressData.table[0].data[i])

    # 预混燃烧过程和扩散燃烧过程
    def Burn(self, Rp, SOC=-5, alpha=1.5, Hu=42700e3, L0=14.3):

        assert 0 <= Rp < 1, "premixed burn fraction value error!!"
        self.Hu = Hu
        self.Rp = Rp
        self.L0 = L0
        self.alpha = alpha
        self.SOC = SOC

        # 循环喷油量
        self.mix.gf = self.mix.M_air(0) / L0 / alpha
        print("Injected fuel per cylcle per cylinder{}mg".format(self.mix.gf * 1.e6))

        self.T2 = self.CompressData.linearInterpolate(SOC, 3)
        self.p2 = self.CompressData.linearInterpolate(SOC, 2)

        # 计算预混燃烧过程
        def fun(x, T):
            print(self.mix.cv(x, T))
            return (Hu - self.mix.u(x, T)) / (alpha * L0 + x) / self.mix.cv(x, T)

        x = 0
        T = self.T2
        step = Rp / 100.
        while x < Rp:
            T += fun(x, T) * step
            x += step
        print("Temperature after premixed burn {}".format(T))

        # 计算扩散燃烧过程
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

        # 计算燃烧终点的曲轴转角
        self.V3 = self.CylGeo.V(0) * T / Ttemp
        self.EOC = self.CylGeo.getFi(self.V3)
        print("End of combustion {}".format(self.EOC))
        self.T3 = T
        # self.p3=self.mix.M_total(1)*self.mix.Rg_gas(1)*self.T3/self.V3
        # print("Pressure after burned {}".format(self.p3))

        return T

    def Expense(self, ke=1.33):
        # 根据燃烧终点的状态去计算膨胀过程缸内温度和压力的变化
        self.ke = ke
        from numpy import arange
        from .Valve import mod
        for i in arange(self.ValveTiming.IVC, self.SOC + 360):
            V = self.CylGeo.V(i)
            p = self.p3 * pow(self.V3 / V, ke)
            T = self.T3 * pow(self.V3 / V, ke - 1)
            self.ExpanseData.append([mod(i), V, p, T, self.mix.M_total(1)])
        self.ExpanseData.doQuickSort(0)

    def pressureReconstructSingeWiebe(self, m):
        from numpy import arange
        from Cylinder import SingleWiebe
        heatRelease = SingleWiebe(self.SOC, self.EOC - self.SOC, m)

        for i in arange(self.SOC, self.EOC):
            p = (1 - heatRelease.X(i)) * self.CompressData.linearInterpolate(i, 2) + heatRelease.X(
                i) * self.ExpanseData.linearInterpolate(i, 2)
            self.Rpressure.append([i, p])

        from numpy import arange
        for i in arange(self.SOC, self.EOC):
            self.data.append([i, self.Rpressure.linearInterpolate(i, 1)])

    def pressureReconstruct(self, heatRelease):
        """
        重构缸压
        :param m: 形状因子
        :return:
        """
        from numpy import arange
        # Sing.plot()
        # 重构缸压
        for i in arange(self.SOC, self.EOC):
            p = (1 - heatRelease.X(i)) * self.CompressData.linearInterpolate(i, 2) + heatRelease.X(
                i) * self.ExpanseData.linearInterpolate(i, 2)
            self.Rpressure.append([i, p])

        from numpy import arange
        for i in arange(self.SOC, self.EOC):
            self.data.append([i, self.Rpressure.linearInterpolate(i, 1)])

    def pressureReconstruct2(self):
        from Cylinder import WibeFunction, SingleWiebe, DoubleWiebe
        from numpy import arange
        comdur = self.EOC - self.SOC
        double = DoubleWiebe(self.Rp, self.SOC, comdur / 2, 3 * comdur / 4)
        double.plot()

        # 重构缸压
        for i in arange(self.SOC, self.EOC):
            # temp = WibeFunction(i, self.SOC, self.EOC - self.SOC, m=m)
            # p = (1 - temp) * self.CompressData.linearInterpolate(i, 2) + temp * self.ExpanseData.linearInterpolate(i, 2)
            p = (1 - double.X(i)) * self.CompressData.linearInterpolate(i, 2) + double.X(
                i) * self.ExpanseData.linearInterpolate(i, 2)
            self.Rpressure.append([i, p])

        from numpy import arange
        for i in arange(self.SOC, self.EOC):
            self.data.append([i, self.Rpressure.linearInterpolate(i, 1)])

    def pit(self, etaTK=0.56, pik=None):
        """
        计算涡前压力
        :param pik: 压气机压比,默认为进气压力除以环境压力
        :param etaTK: 涡轮增压器效率
        :return:
        """
        if pik is None: pik = self.pim / self.p0
        from Compressor import piT
        self.etaTK = etaTK
        from GasProperty import k_exhaust_gu, k_Justi
        # 计算排气门打开时的压力和温度
        # pevo = self.ExpanseData.linearInterpolate(self.ValveTiming.EVO, 2)
        # Tevo = self.ExpanseData.linearInterpolate(self.ValveTiming.EVO, 3)
        pevo = self.ExpanseData.linearInterpolate(180, 2)
        Tevo = self.ExpanseData.linearInterpolate(180, 3)
        print("Temperature at EVO {}".format(Tevo))

        # 等熵膨胀到涡前压力
        # def fun(x):
        #     k = k_exhaust_gu(500)
        #
        #     def Ttfun(kk):
        #         Tt0 = Tevo * pow(x * self.pe0 / pevo, (kk - 1) / kk)
        #         return k_Justi(Tt0, self.alpha) - kk
        #
        #     h = 0.01
        #     while abs(Ttfun(k)) > 1.e-5:
        #         k -= Ttfun(k) / ((Ttfun(k + h) - Ttfun(k - h)) / 2. / h)
        #     # print(k)
        #     self.Tt = Tevo * pow(x * self.pe0 / pevo, (k - 1) / k)
        #     return piT(pik, etaTK, self.Tt, self.T0, self.alpha * 14.3) - x

        # self.pe0

        # 等容膨胀到涡前压力，此时涡轮前温度肯定小于真实的温度，后面来校正
        def fun2(x):
            if self.Tt is None:
                self.Tt = self.pe0 * x / pevo * Tevo
            else:
                pass
            return piT(pik, etaTK, self.Tt, self.T0, self.alpha * 14.3) - x

        x = 2.0
        h = 0.01
        # while abs(fun(x)) > 1.e-5:
        #     x -= fun(x) / ((fun(x + h) - fun(x - h)) / 2. / h)

        while abs(fun2(x)) > 1.e-5:
            x -= fun2(x) / ((fun2(x + h) - fun2(x - h)) / 2. / h)

        print("Pressure ratio of turbine {}".format(x))
        print("pressure before turbine {} bar".format(x * self.pe0 / 1.e5))
        print("Temperature before turbine by isentropic expense {} K".format(self.Tt))

        self.pem = x * self.pe0
        return self.pem

    # def pit2(self, etaTK=0.56, pik=None):
    #     if pik is None: pik = self.pim / self.p0
    #     from Compressor import piT
    #     from GasProperty import k_exhaust_gu, k_Justi

    # 根据涡轮前的状态和涡轮当量流通面积计算流量
    def massFlowRate(self, Area):
        from Cylinder import Volume
        V = Volume(1, self.pem, self.Tt, self.alpha)
        V2 = Volume(1, self.pe0, self.T0, 1.e8)
        from Valve import ValveSimple
        Valve = ValveSimple(Area)
        Valve.connect_to(V, V2)
        return Valve.update(1)

    def gasExchange(self):
        pem = self.pem
        from numpy import arange
        int = self.ValveTiming.int = (self.ValveTiming.EVC + self.ValveTiming.IVC) / 2.
        for i in arange(self.ValveTiming.EVC, int):
            self.data.append([i, self.pim])
        for i in arange(int, self.ValveTiming.IVC):
            temp = xi(i, int, self.ValveTiming.IVC)
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

    def plot(self, original_curve=None):
        self.data.doQuickSort(0)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, figsize=(10, 5))
        ax.plot(self.data.table[0].data, self.data.table[1].data, label="Reconstructed cylinder pressure")
        compre = self.CompressData.slice(self.ValveTiming.IVC, self.ValveTiming.EVO)
        ax.plot(compre.table[0].data, compre.table[2].data, "r-.", label="Compression curve")
        expen = self.ExpanseData.slice(self.SOC, self.ValveTiming.EVO)
        ax.plot(expen.table[0].data, expen.table[2].data, "r-.", label="Expansion curve")

        plt.xlabel(self.data.table[0].ColName + "(" + self.data.table[0].ColUnit + ")")
        plt.ylabel(self.data.table[1].ColName + "(" + self.data.table[1].ColUnit + ")")

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

        ax.scatter(self.ValveTiming.IVC, self.data.linearInterpolate(self.ValveTiming.IVC, 1))
        ax.annotate('IVC %.3g $^\circ$CA' % self.ValveTiming.IVC,
                    xy=(self.ValveTiming.IVC, self.data.linearInterpolate(self.ValveTiming.IVC, 1)), xycoords='data',
                    xytext=(0, 10), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))

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
        if original_curve is not None:
            ax.plot(original_curve.table[0].data, original_curve.table[1].data,
                    label="Original cylinder pressure curve")

        plt.legend()
        plt.tight_layout()
        plt.show()

    def analyze(self, speed=1500, plot=False):
        # 将压缩过程的缸压记录到重构的曲线中
        from numpy import arange
        for i in arange(self.ValveTiming.IVC, self.SOC):
            self.data.append([i, self.CompressData.linearInterpolate(i, 2)])

        # 将燃烧结束到排气门打开这一过程记录
        from numpy import arange
        for i in arange(self.EOC, self.ValveTiming.EVO):
            self.data.append([i, self.ExpanseData.linearInterpolate(i, 2)])

        # 重新排序
        self.data.doQuickSort(0)

        # 循环时间
        tcylce = 30 * 4 / speed

        # 计算PV图
        PV = ArrayTable(3, 0)

        for i in range(self.data.row):
            PV.append(
                [self.data.table[0].data[i], self.CylGeo.V(self.data.table[0].data[i]), self.data.table[1].data[i]])

        # 计算温度
        mass = self.mix.M_total(0)
        from GasProperty import Rg
        temprature = self.CompressData.selectColumns([0, 3]).slice(self.ValveTiming.IVC, self.SOC)
        # temprature.clear()
        #
        # from numpy import arange
        # for i in arange(self.ValveTiming.)

        for i in range(self.Rpressure.row):
            phi = self.Rpressure.table[0].data[i]
            Vtemp = self.CylGeo.V(phi)
            p = self.Rpressure.table[1].data[i]
            temprature.append([phi, p * Vtemp / mass / Rg()])

        ExpanseTemp = self.ExpanseData.slice(self.EOC, self.ValveTiming.EVO)
        for i in range(ExpanseTemp.row):
            temprature.append([ExpanseTemp.table[0].data[i], ExpanseTemp.table[3].data[i]])
        # temprature.plot()

        # PV.plot(2, 1)

        # 单缸指示功
        IndecatedWork = PV.integrate(2, _colx=1)

        # 单缸指示功率
        IndicatedPower = IndecatedWork / tcylce
        print("Indicated power{}".format(IndicatedPower * self.CylGeo.num_of_cylinders))

        # 指示平均有效压力
        IMEP = IndecatedWork / self.CylGeo.displacedVolume()
        print("IMEP={}bar".format(IMEP / 1.e5))

        # 指示热效率
        etait = IndecatedWork / (self.mix.gf * self.Hu)
        print("thermal efficiency {}".format(etait))

        # 迭代求解平均有效压力
        # \[IMEP*{\eta _m}\left( {{C_m},BMEP} \right) = BMEP\]
        BMEP=BMEP2(IMEP,self.CylGeo.bore,self.CylGeo.stroke,speed)
        # def fun(BMEPpara):
        #     return BMEPpara - IMEP * MeEff(self.CylGeo.bore, self.CylGeo.stroke * speed / 30, BMEPpara)

        # BMEP = IMEP;
        # h = IMEP / 1.e2
        # while abs(fun(BMEP) > 1.e-5):
        #     BMEP -= fun(BMEP) / ((fun(BMEP + h) - fun(BMEP - h)) / 2 / h)
        print("BMEP={} bar".format(BMEP / 1.e5))

        etam = MeEff(self.CylGeo.bore, self.CylGeo.stroke * speed / 30, BMEP)
        print("Mechanical efficiency {}".format(etam))

        # 计算有效功率,W
        BrakePower = BMEP * self.CylGeo.totalDisplacedVolume() / tcylce
        print("Brake power {}kW".format(BrakePower / 1.e3))

        # 计算油耗率,kg/h
        B = self.mix.gf * self.CylGeo.num_of_cylinders / tcylce * 3600

        # 计算燃油消耗率，g/(kW*h)，与进气质量和过量空气系数相关
        BSFC = B * 1.e3 / (BrakePower / 1.e3)

        R = thermalUsageIndex(BSFC * etam, self.alpha, self.CylGeo.stroke * speed / 30, self.CylGeo.bore,
                              self.CylGeo.stroke, speed, self.pim, self.Tim)
        self.Tt = exhaustTemperature(BSFC, etam, self.alpha, 1, self.Tim, R)
        self.pit(self.etaTK)
        print("Exhasut temperature by equation {} K".format(self.Tt))

        # 计算单缸流量
        self.massflow = self.mix.M_total(1) / tcylce * self.CylGeo.num_of_cylinders
        print("BSFC={}g/(kW*h)".format(BSFC))

        print("Total intake mass flow rate {} kg/s".format(self.massflow))

        print("maximum pressure {}bar".format(max(self.data.table[1].data) / 1.e5))

        if plot:
            self.plot()

        # print("Total turbine mass flow rate of turbine {}kg/s".format(self.massFlowRate(0.0005)))

        return BMEP, BSFC, self.Tt  # , self.massFlowRate(0.0005)

    def loop(self):
        # 输出平均有效压力，油耗，
        # 气缸几何参数，充量系数，过量空气系数，CA0,燃烧形状因子
        # 调参，排气温度，进气压力，进气温度,涡轮增压器效率
        pass

    def do(self,pim=2.e5,Tim=300,xr=0.1,kc=1.35,phic=0.9,Rp=0.1,alpha=2,SOC=-10,ke=1.37,m=2,etaTK=0.45,speed=2100,plot=False):
        self.compress(pim,Tim,300,xr,kc,phic)
        self.Burn(Rp=Rp,alpha=alpha,SOC=SOC)
        self.Expense(ke)
        self.pressureReconstructSingeWiebe(m)
        self.pit(etaTK)
        self.gasExchange()
        self.analyze(speed,plot=plot)
        return CylinderPressure(self.data,self.CylGeo)


if __name__ == "__main__":

    import json

    CylinderGeometry("WP7").plotVolume()

    from Engine.Valve import ValveSimple

    V1 = Volume(1, 1.e5, 300, 1)
    V2 = Volume(2, 6.e5, 400, 10)
    V3 = Volume(3, 1.e5, 500, 20)
    valve = ValveSimple(1.e-3)
    valve2 = ValveSimple(2.e-3)
    valve3 = ValveSimple(2.e-3)
    valve.connect_to(V1, V2)
    valve2.connect_to(V2, V3)
    valve3.connect_to(V3, V1)

    dt = 0.01
    t = 0
    while t < 10:
        valve.update(t)
        valve2.update(t)
        valve3.update(t)

        V1.m += V1.dm() * dt
        V1.T += V1.dT() * dt
        V1.alphak += V1.dAlpha() * dt

        V2.m += V2.dm() * dt
        V2.T += V2.dT() * dt
        V2.alphak += V2.dAlpha() * dt

        V3.m += V3.dm() * dt
        V3.T += V3.dT() * dt
        V3.alphak += V3.dAlpha() * dt
        # V3.m += V3.dm() * dt
        # V3.p += V3.dp() * dt

        V1.record(t)
        V2.record(t)
        V3.record(t)
        # V3.record(t)
        print("mass={}".format(V1.m + V2.m))

        t += dt
        print(V1.m)
        print(V1.p)
        # os.system("pause")

    valve.data.plot(2)
    V2.data.plot(4)
    V1.data.plot(4)
    V2.data.plot(3)
    V1.data.plot(3)
    V1.data.plot(2)
    V2.data.plot(2)
    # V2.data.plot(3)
    # V1.data.plot(3)

    # V3.data.plot(2)
    # print(V1.data.table[2].data[-1],V2.data.table[2].data[-1],V3.data.table[2].data[-1])
