
import sys

sys.path.append("../Properties of working fluids/")
sys.path.append("../Data/")
sys.path.append("../Pipe/")
sys.path.append("../")
from Cylinder import CylinderGeometry


if __name__=="__main__":
    import json
    with open("IdealMillerCycle.json", mode='r', encoding='UTF-8') as f:
        setting = json.load(f)

    if setting["engine type"] is None:
        cyl=CylinderGeometry(setting["bore,mm"],setting["stroke,mm"],setting["connecting rod length,mm"],setting["compression ratio"],setting["number of cylinders"])

    else:
        from pandas import read_excel
        data=read_excel(setting["data file"],index_col="机型").loc[setting["engine type"]]
        cyl=CylinderGeometry(data["bore,mm"]*1.e-3,data["stroke,mm"]*1.e-3,data["connecting rod length,mm"]*1.e-3,data["compression ratio"],data["number of cylinders"])

    cyl.plotVolume()