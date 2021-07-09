import os
import sys

sys.path.append("../Turbocharging")
sys.path.append("../Valve")
from Cylinder.Cylinder import *
from Valve.valve import *
from pandas import read_excel

if __name__ == "__main__":



    # 气缸容积
    data = read_excel("../Data/机型数据.xlsx", index_col="机型").loc["WP7"]
    cyl = CylinderGeometry("WP7")

    # 放热率
    hrr = DoubleWiebe(0.21, -5, 50, 50)
    # hrr.plot()

    # 传热
    ht = HeatTransfer(500,500,500)

    n = 3000

    C = cylinder(-180, 3.e5, 300, 100, n, cyl, hrr, ht, Inj(230e-6))

    Vol=Volume(1,1.e5,500,10)
    Vol2 = Volume(1, 2.5e5, 300, 1.e5)
    # from valve import Valve
    valvelift = ArrayTable()
    valvelift.readCSVFile("../Valve/IntakeValveLift.csv")
    val2 = Valve(40e-3, 120, 400, valvelift, 0,1.1)
    val1=Valve(50e-3, -400, -150, valvelift, 0)
    val2.connect_to(C,Vol)
    val1.connect_to(Vol2,C)

    # val1.plotValveLift(val2)
    # val2.valve_lift.plot()
    cycle=10

    import time as t

    now = t.time()

    dFi=2
    while C.Fi<-180+720*cycle:
        try:
            C.Fi+=dFi
            # print(val2.massFlowRate(C.Fi,C.N))
            val2.massFlowRate(C.Fi, C.N)
            val1.massFlowRate(C.Fi,C.N)
            C.m+=C.dm()*dFi
            C.T+=C.dT()*dFi
            C.alphak+=C.dAlpha()*dFi
            # C.update()
            C.RecordThisStep()
            Vol.record(C.Fi)
        except:
            C.plot()
    # C.data.openWithProgram()
    print("calculate time {}".format(t.time()-now))
    print("real time {}".format(120/n*cycle))
    C.plot()
    # Vol.data.plot(2)
