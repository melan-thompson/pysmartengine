from Cylinder import *
from ArrayTable import *

Table=ArrayTable()
Table.readCSVFile("wp7缸压0.15.csv")

P=CylinderPressure(Table)

P.FFTTrasform(2200).plot()