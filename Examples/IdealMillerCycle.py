
import sys

sys.path.append("../Properties of working fluids/")
sys.path.append("../Data/")
sys.path.append("../Pipe/")
sys.path.append("../")
sys.path.append("../Turbochaging/")
from Cylinder.Cylinder import CylinderGeometry
from Cylinder.Cylinder import MillerCycle


if __name__=="__main__":
    import json
    with open("IdealMillerCycle.json", mode='r', encoding='UTF-8') as f:
        setting = json.load(f)

    if setting["engine type"] is None:
        cyl=CylinderGeometry(setting["engine type"])

    else:
        from pandas import read_excel
        data=read_excel(setting["data file"],index_col="机型").loc[setting["engine type"]]
        cyl=CylinderGeometry(setting["engine type"])

    cylcle=MillerCycle(cyl,setting["environment pressure,Pa"],setting["environment temperature,K"])
    cylcle.initCompressor(setting["compressor efficiency"],setting["compressor pressure ratio"])
    cylcle.intake(setting["intake valve close timing,bBTC"])
    cylcle.adiabaticompress()
    cylcle.preBurn(setting["premixed burned fuel ratio"],setting["excess air fuel ratio"])
    cylcle.disffBurn()
    cylcle.expansion()
    cylcle.supercriticalExhaust(setting["turbo efficiency"])
    cylcle.subcritical()
    cylcle.analyze(plot=False)
    cylcle.plot(1)