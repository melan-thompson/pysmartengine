# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

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
#from GasProperty import *
from Compressor import *

from FileManipulate import *
from Pipe import QPM
from Algorithm import GP
from MEEM import *
from Valve import *

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # T=ArrayTable(1,0)
    # T.append([6])
    # T.append([1])
    # T.append([2])
    # T.append([7])
    # T.append([9])
    # T.append([3])
    # T.append([4])
    # T.append([5])
    # T.append([10])
    # T.append([8])
    # T.doQuickSort(0)
    # T.openWithProgram()

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

    C = CylinderGeometry(100e-3, 100e-3, 152e-3, 18)
    print(massAir(C.displacedVolume(),2.7649e5,318))
    I=IdealCylcle(C)
    I.compress(-110,xr=0.07,Tim=318,Tr=1000,pim=2.7649e5)
    # I.CompressData.plot(2)
    I.Burn(Rp=0.01,alpha=1.7,SOC=1.5)
    I.Expense()
    I.pressureReconstruct()
    I.gasExchange(pem=3.33e5)
    I.analyze()
    I.plot()
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

    def fun(t,y):
        from math import exp
        return y*exp(t)

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
    # cylce=MillerCycle(C)
    # cylce.initCompressor(0.9,2)
    # cylce.intake(100)
    # cylce.adiabaticompress()
    # cylce.preBurn(0.1)
    # cylce.disffBurn()
    #
    # # cylce.premixBurn(1)
    # # cylce.DisffusionBurn()
    # cylce.expansion()
    # cylce.supercriticalExhaust(0.8)
    # cylce.subcritical()
    # # cylce.data.openWithProgram()
    # cylce.data.animation(2,1)
    # cylce.data.plot(3,1)


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
