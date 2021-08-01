import os
import sys

sys.path.append("../..")
from FluidProperties.GasProperty import *
from Cylinder.Cylinder import MeEff, thermalUsageIndex, exhaustTemperature

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

        from ArrayTable import ArrayTable
        Table = ArrayTable()
        Table.readCSVFile("EnergyBalance.csv")
        qw = Table.GPR([bore * 1.e3, n, self.pme / 1.e5], [1, 4, 5], 9, trained=False)[0][0] / 1.e2
        print(qw)

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
            os.system("start test2.html")
        # snakey = Sankey()
        # snakey.add(
        #     "",
        #     self.nodes,
        #     self.links
        # )
        # snakey.render(path="chart.pdf", delay=3)


if __name__ == "__main__":
    E = EnergyBalance(160e-3, 160e-3, 1500, 216.4, 1105.76, 75 + 273.15, 1.5e5, 3, 12,1.1,2,0.5,20+273.15)
    E.plot()
