from GasProperty import *

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


def piK(Tt, pit, etat, etak, T0=273.15 + 20, air_fuel_ratio=23.85, ma=1, me=1, etam=1):
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
    temp = (cp_Justi(T0) * T0*air_fuel_ratio) / ((air_fuel_ratio + 1) * cp_exhaust_gu(Tt) * Tt * etatk)
    temp2 = pow(pik, (gammak - 1) / gammak) - 1
    return pow(1 - temp * temp2, -gammat / (gammat - 1))

# class Compressor:
#     def __init__(self, _map):
#         self.map = _map
