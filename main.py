# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import math
import sys

sys.path.append("./Properties of working fluids/")
sys.path.append("./Data/")
sys.path.append("./Pipe/")
sys.path.append("./Algorithm/")
sys.path.append("./Mean value engine models/")
from ArrayTable import *


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


from Cylinder import *
# from GasProperty import *
from Compressor import *

from FileManipulate import *
from Pipe import QPM
from Algorithm import GP
from MEEM import *
from Valve import *
from GP import *

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # arraytable=ArrayTable()
    # arraytable.readExcelFile("油耗.xlsx",sheetname=0,_header=[0,1])
    # arraytable.groupPlot()
    # arraytable.surfaceFit()
    # arraytable.scatter3D()
    # arraytable.openWithProgram()

    # import pandas as pd
    # import matplotlib.pyplot as plt
    #
    # plt.rcParams['font.sans-serif'] = ['SimHei']

    import matplotlib as mpl
    import numpy as np


    # Make a modified version of cdict with some transparency
    # in the middle of the range.cdict1['alpha'] = ((0.0, 1.0, 1.0),
    #   (0.25,1.0, 1.0),
    # (0.5, 0.3, 0.3),
    #   (0.75,1.0, 1.0),
    # (1.0, 1.0, 1.0))


    # pp=pd.DataFrame([[1,2],[4,5]],columns=['a','b'])
    # pp.plot.scatter(x='a',y='b')
    # plt.show()
    # sheet=pd.read_excel("功率计算.xlsx",sheet_name=None)
    # print(list(sheet.keys()))

    # df = pd.read_excel("机型数据.xlsx", sheet_name=0, header=0, usecols=range(85))
    # result = df.columns.values
    # print(result)
    # group=df.groupby("应用")
    # cleaned=df.dropna()
    # name=list()
    # for each,ii in group:
    #     name.append(each)

    data = pd.read_excel("油耗.xlsx")


    # name = list()
    # group = data.groupby("转速")
    # for each, ii in group:
    #     name.append(each)
    # ii = 0
    # ax=None
    # for each in name:
    #     ax = group.get_group(each).plot(x="扭矩", y="燃油消耗率", label=each, ax=ax)
    #     ii += 1
    # plt.show()
    #
    # plt.figure()
    # ax = None
    # df.plot.scatter(x="标定转速，rpm", y="平均有效压力,Mpa")
    # ii=0
    # for each in name:
    #     ax=group.get_group(each).plot.scatter(x="缸径,mm",y="冲程,mm",label=each,ax=ax,c=cnames[list(cnames.keys())[ii]])
    #     ii+=1
    # plt.scatter(30,50,color='r')
    # plt.tight_layout()
    # # plt.legend()
    # plt.show()
    # df.plot.show()
    # print(result[3], result[4])
    # df.to_excel("you.xlsx",sheet_name="how")
    # from os import system
    # system("start you.xlsx")
    # system("pause")
    # print(df.columns.values)
    # result=df.columns.values
    # print(type(list(df[result[2]])))
    # print(df.head())
    # print(TurbochargePressure(pme=14.7891e5,T_im=335.71,eta_et=0.3695,phi_a=1.668,VE=0.8637))
    # Gauss2DSample()

    # def abs(x):
    #     return x
    #
    # import torch as t
    # from torch.autograd import Variable as V
    # x=V(t.ones(1),requires_grad=True)
    # y=abs(x)
    # y.backward()
    # print(x.grad)

    # print(exhaustTemperature())

    # def combustionEfficiency(alpha):
    #     return 0.98 * min(1, 1.2 * alpha - 0.2)
    #
    #
    # table = ArrayTable(2,0)
    # for i in range(0,7):
    #     table.append([i,combustionEfficiency(i)])
    # table.plot()
    # table.plotfun(combustionEfficiency, 1, 7)

    valveTiming = ValveDesign()
    valveTiming.changeTiming(IVO=-377, IVC=-154, EVO=125, EVC=375)
    # valveTiming.plot()

    WP7 = CylinderGeometry(108e-3, 130e-3, 209.7e-3, 18, 6)



    pressure = ArrayTable()
    pressure.readCSVFile("wp7缸压0.15.csv")

    pre = CylinderPressure(pressure)
    # pre.plot(valveTiming)
    # pre.slice(-60,60)
    # pre.polyTropicIndex(WP7,plot=True)
    # pre.netHeatReleaseRate(WP7, -110, plot=True)
    # pre.slice(-110,-6)
    # pre.LogP_LogV(WP7).polyFit().plot([1,2])
    # pre.LogP_LogV(WP7, plot=True, ivc=-154, evo=125)
    # pre.data.polyFit(_deg=5).plot([1,2])
    # pre.EquivalentIsotropicIndex(WP7,-90).slice(-60,-20).plot()
    # pre.PVDiagram(WP7).Log_plot()
    # pre.startOfCombustion(2)
    # pre.plot()

    # C=ValveDesign(TDC=0)
    # C.plot()
    # from numpy import linspace
    # def divideCV(start=0,end=10,numberofCV=10):
    #     from numpy import linspace
    #     boundary=linspace(start,end,numberofCV+1)
    #     centerofCV=[(boundary[i]+boundary[i+1])/2. for i in range(len(boundary)-1)]
    #     return centerofCV,boundary
    #
    # node,B=divideCV(0,10,10)
    # print(node)
    # print(B)
    from pandas import read_excel

    # data=read_excel("input.xlsx",header=None,sheet_name=0)
    # print(data[1])

    # BSFCexample(2000).plot()
    # TransposeCSVFile("ENGINEDATA.csv")

    # import sqlite3
    #
    # conn = sqlite3.connect("temp.db")
    # cursor = conn.cursor()
    # cursor.execute("INSERT INTO temptable VALUES (1,2,3,4,5),(5,6,89,12,45);")
    # conn.commit()
    # conn.close()
    # table=ArrayTable(2,0)
    # from numpy import arange
    # for i in arange(0.5,10.1):
    #     table.append([i,IdealMillerCylcle(Rp=0.06,alpha=2,eps=16,epsm=14,etaTK=0.5,pik=i,T2=350)])
    # table.plot()
    # print(nozzle(1.e5))
    # table.plotfun(x,-120,60).plot()
    # table.plotfun(nozzle, 1.e5,1.e6).plot(0, 1)

    ###################缸压重构######################
    C = CylinderGeometry(100e-3, 100e-3, 152e-3, 18)
    print(massAir(C.displacedVolume(), 2.7649e5, 318))
    I = IdealCylcle(WP7, valveTiming)
    I.compress(xr=0.07, Tim=318, Tr=1000, pim=200.9e3, kc=1.38, phic=1)
    # I.CompressData.plot(2)
    I.Burn(Rp=0.01, alpha=1.63, SOC=-6)
    I.Expense(ke=1.29)
    I.pressureReconstruct(m=1.9)
    I.pit(p0e=1.e5 + 0.15e5, etaTK=0.6)
    I.gasExchange()
    I.data.doQuickSort()
    I.analyze(plot=True)
    # I.plot()
    # I.data.openWithProgram()
    pressure.setUnitToSI()
    # pressure.plot()

    # I.ExpanseData.plot(3)
    I.data.compareResultPlot([pressure])


    # I.analyze()
    # I.plot()

    # I.Rpressure.plot()
    # I.ExpanseData.plot(2)
    # I.ExpanseData.compareResultPlot([I.CompressData,I.GasExchange])
    # I.GasExchange.compareResultPlot([I.CompressData])
    # I.CompressData.plot(2)
    # table.plotfun(openEnd2).plot()
    # table.readCSVFile("机型数据库.csv",1,typename="string")
    # table.translate()
    # table.openWithProgram()
    # table.insertIntoSQLDB("EngineDB.db","CylinderGeometry")
    # table.openWithProgram()
    # table.readSQLiteTable("EngineDB.db","EngineBasicPara")
    # table.readSQLiteTable("EngineDB.db","CylinderGeometry")
    # table.selectValidData(isNumber=False).openWithProgram()
    # table.selectValidData([1,10,12,13]).selectValidData([0,1,2,3],isNumber=False).openWithProgram()
    # table.openWithProgram()
    # table.colorScatter(11,13,14)
    # table.readCSVFile("统计值.csv")
    # table.setUnitToSI()
    # table.writeToSQLite("FuelConsumptionResult.db","FuelConsumption")

    # table.colorScatter()
    # table.interpolate2D()
    # table.scatter(1)

    def fun(t, y):
        from math import exp
        return y * exp(t)


    def fun1(x):
        import math
        return math.sin(x)
    # table=ArrayTable()
    # table.plotfun(fun1,-10,10).plot()
    # solver=ODEsolver(fun,0,1)
    # solver.implicitEuler(2,0.01).compareResultPlot([solver.Rungekutta4(2),solver.forecastCorrection(2),solver.eulersolver(2)])
    # solver.Rungekutta4(2).plot()
    # solver.forecastCorrection(2).plot()
    # solver.eulersolver(2).plot()

    # import numpy as np
    # print((np.expand_dims(np.array([1,2,3]),1)-np.expand_dims(np.array([1,2,3]),0)**2))
    # polt2DGaussianDistribution([-3,3])
    # print(piT(2.875456,0.5987,884.017,295.56,23.85))
    # print(piT(2.5794, 0.5924, 872.32, 296.40, 23.63))
    # print(piT(2.23094, 0.57495, 866.023, 297.04, 22.94))
    # table=ArrayTable(3,0)
    # table.setTableHeader(["piK","piK","piT"])
    # table.setTableUnit(["/","/","/"])
    # from numpy import arange
    # for i in arange(1.,2,0.01):
    #     table.append([i,i,piT(i,0.6,550)])
    # table.plot([1,2])

    # C=CylinderGeometry(100e-3,100e-3,160e-3,16)
    # # C.plotVolume()
    cylce=MillerCycle(C)
    cylce.initCompressor(0.9,2)
    cylce.intake(100)
    cylce.adiabaticompress()
    cylce.preBurn(0.1)
    cylce.disffBurn()

    cylce.premixBurn(1)
    cylce.DisffusionBurn()
    cylce.expansion()
    cylce.supercriticalExhaust(0.8)
    cylce.subcritical()
    # cylce.data.openWithProgram()
    cylce.data.animation(2,1)
    cylce.data.plot(3,1)

    # tt=1
    # import numpy as np
    # import matplotlib.pyplot as plt
    # from matplotlib.animation import FuncAnimation
    #
    # fig, ax = plt.subplots()
    # #plt.autoscale(enable=True,axis="both",tight=True)
    # #xdata, ydata = [], []
    # ln, = ax.plot([], [], 'r-', animated=False)
    #
    #
    #
    # def init():
    #     ax.set_xlim(-1000, 1000)
    #     ax.set_ylim(0, 2000)
    #     return ln,
    #
    #
    # def update(t):
    #     solution=AnalyticalSolution(t, [1, 1200, 0.5e6], [0, 300, 0.1e6],[-1000,1000])
    #     #ax.set_xlim(min(solution.table[0].data), max(solution.table[0].data))
    #     #ax.set_ylim(0.9 * min(solution.table[4].data), 1.1 * max(solution.table[4].data))
    #     ln.set_data(solution.table[0].data, solution.table[4].data)
    #
    #     return ln,
    #
    #
    # ani = FuncAnimation(fig, update, frames=np.linspace(1.e-5, 1, 128),
    #                     init_func=init, blit=True)
    # plt.show()
    # shoc=ShockTube([1, 1200, 0.5e6], [0, 300, 0.1e6])
    # shoc.animation(3)
    # C=CylinderGeometry(108.e-3,130.e-3,209.7e-3,16)
    # temp=ArrayTable()
    #
    # temp.readCSVFile("Input-4.csv","示功图")
    # temp.setTableHeader(["Crank angle", "Cylinder pressure"])
    # temp.setTableUnit(["°CA", "bar"])

    # temp.diff(1)
    # temp.diff(2)
    # temp.diff(3)
    # temp.plot([3,4])
    # CylPre=CylinderPressure(temp)
    # CylPre.smooth("five points three times smooth").openWithProgram()
    # CylPre.data.plot()
    # CylPre.data.plot([1,2])
    # CylPre.startOfCombustion(2)
    # CylPre.slice(-50,50)
    # CylPre.data.diff(1)
    # CylPre.data.diff(2)
    # CylPre.data.diff(3)
    # CylPre.data.plot(3)
    # CylPre.data.plot(4)
    # CylPre.Temperature(C).plot()
    # CylPre.PVDiagram(C).plot()
    # temp.plot()

    # p=CylinderPressure(temp)
    # p.slice(-100,100)
    # # p.data.plot()
    # p.FFTTrasform(1500).plot()
    #
    # temp.diff(1)
    # temp.diff(2)
    # temp.diff(3)
    # temp.Log_plot(1)
    #
    # Geo=CylinderGeometry(108e-3,130e-3,209e-3,16)
    # # Geo.plotVolume()
    # temp=ArrayTable()
    # temp.readCSVFile("Input-4.csv","示功图")
    # pre=CylinderPressure(temp)
    # pre.sliceToComAndExhaust()
    # pre.ployTropicIndex(Geo).plot()
    #
    # pre.ployTropicIndex(Geo).openWithProgram()
    # temp.diff(1)
    # temp.animation(2)
    # def fun(T):
    #     return 1.338-6.0e-5*T+1.0e-8*T*T
    # temp=ArrayTable(2,0)
    # temp.setTableHeader(["Temperature","$\gamma$"])
    # temp.setTableUnit(["K","/"])
    # for i in range(300,3000):
    #     temp.append([i,fun(i)])
    # temp.plot()
    #
    # laxWendroff1step(1.e-4,[1, 1200, 0.5e6], [0, 300, 0.1e6]).writeToSQLite()
    #
    # a=Node()
    # a.init(u=1,p=5.e5,T=1200)
    # print(a.U)
    # print(a.F())
    # print(a.e)
    # a.solve()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
