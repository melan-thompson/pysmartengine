# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

from Cylinder import *
#from GasProperty import *
from Compressor import *
from ShockTube import *
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
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



    C=CylinderGeometry(100e-3,100e-3,160e-3,16)
    # C.plotVolume()
    cylce=MillerCycle(C)
    cylce.initCompressor(0.9,2)
    cylce.intake(100)
    cylce.adiabaticompress()
    cylce.preBurn(0.1)
    cylce.disffBurn()

    # cylce.premixBurn(1)
    # cylce.DisffusionBurn()
    cylce.expansion()
    cylce.supercriticalExhaust(0.8)
    cylce.subcritical()
    # cylce.data.openWithProgram()
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
    # temp=ArrayTable()
    # temp.readCSVFile("Input-4.csv","示功图")
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




# See PyCharm help at https://www.jetbrains.com/help/pycharm/
