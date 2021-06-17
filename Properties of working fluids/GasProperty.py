def u_Justi(T, AFAK=1.e8, Tref=273.15) -> float:
    """
    Justi公式
    :param T:Temperature of gas,废气温度(K)
    :param AFAK:Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :param Tref: Reference temperature，参考温度(K)
    :return:比热力学能(J/kg)
    """
    temp1 = 489.6 + 46.4 / pow(AFAK, 0.93)
    temp2 = 7.768 + 3.36 / pow(AFAK, 0.8)
    temp3 = 0.0975 + 0.0485 / pow(AFAK, 0.75)
    temp4 = 1356.8 + temp1 * (T - Tref) * 1.e-2 + temp2 * pow(T - Tref, 2) * 1.e-4 - temp3 * pow(T - Tref, 3) * 1.e-6
    return 0.1445 * temp4 * 1.e3  # J / kg


def cv_Justi(T, AFAK=1.e8, Tref=273.15):
    """
    由Justi公式算出的比定压热容
    :param T: Temperature of gas,废气温度(K)
    :param AFAK: Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :param Tref: Reference temperature，参考温度(K)
    :return: 定压比热容(J/(kg*K))
    """
    h = 1.e-3
    return (u_Justi(T + h, AFAK, Tref) - u_Justi(T - h, AFAK, Tref)) / 2 / h


def cv_mean_Justi(start, end, AFAK=1.e8, step=0.1):
    if start > end:
        raise Exception("T_begin must be less than T_end when calculating mean cv!!!")
    result = 0
    import numpy as np
    for T in np.arange(start, end, step):
        result += step / 2. * (cv_Justi(T, AFAK) + cv_Justi(T + step, AFAK))
    return result / (end - start)


def cp_mean_Justi(start, end, AFAK=1.e8, step=0.1, ):
    if start > end:
        raise Exception("T_begin must be less than T_end when calculating mean cp!!!")
    result = 0
    import numpy as np
    for T in np.arange(start, end, step):
        result += step / 2. * (cp_Justi(T, AFAK) + cp_Justi(T + step, AFAK))
    return result / (end - start)


def k_mean_Justi(start, end, AFAK=1.e8, step=0.1):
    return cp_mean_Justi(start, end, AFAK, step) / cv_mean_Justi(start, end, AFAK, step)


def Rg(AFAK=1.e8):
    """
    J.keenan & J.Kaye formula, used for ideal gas
    将废气看作理想气体时计算的，气体常数，与温度无关
    :param AFAK: Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :return: 气体常数，与温度无关// J / (kg * K)
    """
    return 9.81 * (29.2647 - 0.0402 / AFAK)


def h_Justi(T, AFAK=1.e8, Tref=273.15):
    """
    Justi公式和J.keenan & J.Kaye公式算出的比焓，适用于理想气体
    :param T:Temperature of gas,废气温度(K)
    :param AFAK:Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :param Tref: Reference temperature，参考温度(K)
    :return:比焓(J/kg)
    """
    return u_Justi(T, AFAK, Tref) + Rg(AFAK) * T


def k_Justi(T, AFAK=1.e8):
    """
    Justi公式和J.keenan & J.Kaye公式算出的比热容比
    :param T:Temperature of gas,废气温度(K)
    :param AFAK:Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :return:比热容比
    """
    return 1 + Rg(AFAK) / cv_Justi(T, AFAK);


def cp_Justi(T, AFAK=1.e8, Tref=273.15):
    """
    Justi公式和J.keenan & J.Kaye公式算出的定压比热容
    :param T:Temperature of gas,废气温度(K)
    :param AFAK:Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :return:定压比热容((J/(kg*K)))
    """
    return cv_Justi(T, AFAK, Tref) * k_Justi(T, AFAK)


def M(AFAK=1.e8):
    """
    J.keenan & J.Kaye公式算出的相对分子质量
    :param AFAK:Generalized excess air fuel ratio,广义过量空气系数，1.e8时为空气
    :return:相对分子质量g / mol
    """
    return 28.9705 + 0.0403 / AFAK


def cv_air_gu(T):
    """
    顾闳中二次三项式定压
    :param T:Temperature of air, 空气温度(K)
    :return:空气比定压比热容(J/(kg·K))
    """
    Ma = 28.96  # 空气的相对分子质量
    return (4.678028662 + 6.872284154e-4 * T - 6.06834439e-8 * pow(T, 2)) * pow(10, 3) * 4.187 / Ma


def cv_exhaust_gu(T, alpha=1):
    """
    顾闳中教授废气定压比热公式
    :param T: Temperature of air, 空气温度(K)
    :param alpha: 空气比定压比热容(J/(kg·K))
    :return:
    """
    Me = M(alpha)
    return (19.893 + 5.020e-3 * T - 5.959e-7 * T ** 2) / Me * 1.e3
    # return (4.751276526 + 1.19900582e-3 * T - -1.42321698e-7 * pow(T, 2)) * pow(10, 3) * 4.187 / Me


def mucps(Ts):
    # \left( {\mu {c_p}} \right)_s = 27.59 + 0.0025{T_s}
    """
    计算排温时用到的空气平均定压摩尔比热容
    :param Ts:进气管中的温度(K)
    :return:空气平均定压摩尔比热容,kJ/(kmol * K)或者J/(mol * K)
    """
    return 27.59 + 2.5e-3 * Ts


def mucpT(Tt, alpha=1.1, phis=1):
    # {\left( {\mu {c_p}} \right)_T} = 8.315 + \frac{{20.47 + \left( {\alpha {\varphi _s} - 1} \right) \times
    # 19.26}}{{\alpha {\varphi _s}}} + \frac{{3.6 + \left( {\alpha {\varphi _s} - 1} \right) \times 2.51}}{{\alpha {
    # \varphi _s} \times {{10}^3}}}{T_T}
    """
    计算排气温度时，温度为$T_T$时废气平均定压摩尔比热容$[kJ/(kmol K)]$
    :param Tt:排气温度(K)
    :param alpha:过量空气系数
    :param phis:扫气系数
    :return:废气平均定压摩尔比热容$[kJ/(kmol K)]$
    """
    return 8.315 + (20.47 + (alpha * phis - 1) * 19.26) / (alpha * phis) + (3.6 + (alpha * phis - 1) * 2.51) / (
            alpha * phis * 1e3) * Tt


def k_exhaust_gu(T, alpha=1):
    return Rg(alpha) / cv_exhaust_gu(T, alpha) + 1


def cp_exhaust_gu(T, alpha=1):
    return cv_exhaust_gu(T, alpha) + Rg(alpha)


def m_fuel(m_total, alphak, L0=14.3):
    # {g_f}{x_k} = \frac{{{G_s}}}{{{L_0}{\alpha _k} + 1}}
    """
    由缸内总质量和广义过量空气系数计算燃油量
    :param m_total: 缸内总质量(kg)
    :param alphak: 广义过量空气系数
    :param L0: 燃空当量比
    :return: 形成混合所需要的燃油量
    """
    if alphak < 1:
        raise Exception("excess air fuel ratio can not less than 1")
    return m_total / (L0 * alphak + 1)


class DieselMixture:
    def __init__(self, L0=14.3):
        self.L0 = L0
        # self.gf = gf
        from ArrayTable import ArrayTable
        self.data = ArrayTable(8, 0)
        self.data.setTableHeader(["已燃百分比",
                                  "缸内气体总质量", "空气质量", "废气质量", "已喷燃油量",
                                  "广义过量空气系数", "残余废气系数",
                                  "气体常数"])
        self.data.setTableUnit(["Fraction",
                                "kg", "kg", "kg", "kg",
                                "", "",
                                "J/(kg·k)"])

    def init_With_Ma_r(self, ma, r, gf):
        self.ma = ma
        self.r = r
        self.gf = gf
        if ma < self.L0 * self.gf:
            print("Charged air is too few to burn all the fuel!!")
            raise Exception("There are must at least {}(kg) air but {} is provided".format(self.L0 * self.gf, ma))
        if r < 0 or r > 1:
            raise Exception("residual gas fraction can not less than 0 and greater than 1!!")
        self.xr = ma * r / (1 - r) / (1 + self.L0) / self.gf
        self.__initPropertyTable()

    def init_With_Mc_r(self, mc, r, gf):
        self.init_With_Ma_r((1 - r) * mc, r, gf)

    def init_With_Mc_AFA(self, mc, afa, gf):
        self.gf = gf
        if afa < 1:
            raise Exception("Excess air fuel ratio must be greater than or equal to 1")
        mr = (1 + self.L0) / (self.L0 * afa + 1) * mc
        self.ma = self.L0 * (afa - 1) / (self.L0 * afa + 1) * mc
        if self.ma < self.L0 * gf:
            print("Charged air is too few to burn all the fuel!!")
            raise Exception("There are must at least {}(kg) air but {} is provided".format(self.L0 * self.gf, self.ma))
        self.r = mr / mc
        self.xr = self.ma * self.r / (1 - self.r) / (1 + self.L0) / self.gf
        self.__initPropertyTable()

    def AFAK(self, x):
        if x < 0: x = 0
        if (self.xr + x) < 1.e-8: return 1.e8
        return (self.L0 * self.gf * self.xr + self.ma) / (self.L0 * self.gf * (self.xr + x))

    def M_total(self, x):
        if x < 0: x = 0
        if (self.xr + x) < 1.e-8:
            return self.ma
        return self.gf * (self.xr + x) * (1 + self.AFAK(x) * self.L0)

    def M_exhaust(self, x):
        return self.M_total(x) * (1 + self.L0) / (self.L0 * self.AFAK(x) + 1)

    def M_air(self, x):
        return self.M_total(x) - self.M_exhaust(x)

    def Xk(self, x):
        return self.M_total(x) / (1 + self.L0 * self.AFAK(x)) / self.gf

    def cv(self, x, T):
        return cv_Justi(T, self.AFAK(x))

    def cp(self, x, T):
        return cp_Justi(T, self.AFAK(x))

    def u(self, x, T):
        return u_Justi(T, self.AFAK(x))

    def h(self, x, T):
        return h_Justi(T, self.AFAK(x))

    def U(self, x, T):
        return self.M_total(x) * u_Justi(T, self.AFAK(x))

    def Rg_gas(self, x):
        return Rg(self.AFAK(x))

    def k(self, x, T):
        return k_Justi(T, self.AFAK(x))

    def __initPropertyTable(self):
        from numpy import arange
        for i in arange(0, 1, 0.01):
            self.data.append([i,
                              self.M_total(i), self.M_air(i), self.M_exhaust(i), self.gf * i,
                              self.AFAK(i), self.M_exhaust(i) / self.M_total(i),
                              self.Rg_gas(i)])


def cp_R(species="H2O", T=3000):
    from pandas import read_excel
    from math import pow
    data = read_excel("Coefficient_for_species_1000_5000.xlsx", index_col="Species")
    result = 0
    coeff = data.loc[species]
    print(coeff)
    for i in range(5):
        result += coeff[i] * pow(T / 1000, i)
        # print(coeff[i])
    # result+=coeff[3]
    # result+=coeff[4]*(T/1000)**5
    return result


if __name__ == "__main__":
    print(cp_R())
