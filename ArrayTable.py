import sys

sys.setrecursionlimit(100000)


class PhsicalVarList:
    def __init__(self, data=list(), ColName="EmptyName", ColUnit="EmptyUnit"):
        self.ColName = ColName
        self.ColUnit = ColUnit
        self.data = data

    def show(self):
        print(self.ColName)
        print(self.ColUnit)
        for each in self.data:
            print(each)

    def plot(self):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.plot(self.data)
        plt.xlabel("Index of data")
        plt.ylabel(self.ColName + "(" + self.ColUnit + ")")

        plt.grid()
        plt.show()

    def findMaxValue(self):
        return max(self.data)

    def findMinValue(self):
        return min(self.data)

    def scaleToRange(self, r=[0, 1]):
        minvalue = self.findMinValue()
        maxvalue = self.findMaxValue()
        result = [0] * len(self.data)
        for i in range(len(self.data)):
            result[i] = r[0] + (self.data[i] - minvalue) * (r[1] - r[0]) / (maxvalue - minvalue)
        return result

    def mean(self):
        return sum(self.data) / len(self.data)

    def Var(self):
        meanvalue = self.mean()
        temp = [(each - meanvalue) ** 2 for each in self.data]
        return sum(temp) / (len(self.data) - 1)

    def moveToCenter(self):
        meanvalue = self.mean()
        result = [self.data[i] - meanvalue for i in range(len(self.data))]
        self.data = result

    def __mul__(self, num):
        self.data = [each * num for each in self.data]

    def __add__(self, num):
        self.data = [each + num for each in self.data]


class ArrayTable:
    def __init__(self, col=0, row=0):
        self.col = col
        self.row = row
        self.table = list([PhsicalVarList([None for j in range(row)]) for i in range(col)])

    def show(self):
        for i in range(len(self.table)):
            print("{:<20}".format(self.table[i].ColName), end="")
        print("\n")
        for i in range(len(self.table)):
            print("{:<20}".format(self.table[i].ColUnit), end="")
        print("\n")
        for i in range(self.row):
            for j in range(len(self.table)):
                print("{:<20}".format(self.table[j].data[i]), end="")
            print("\n")

    def append(self, listToAppend):
        if len(listToAppend) != self.col:
            raise Exception(
                "List to append does not have the same lenth as the table, table has {0} column, but the list has {1} "
                "columns".format(
                    self.col, len(listToAppend)))
        else:
            for i in range(self.col):
                self.table[i].data.append(listToAppend[i])
            self.row += 1

    def appendColumn(self, _data):
        if type(_data) == list:
            for each in _data:
                self.table.append(each)
            self.col += len(_data)
        else:
            self.table.append(_data)
            self.col += 1

    def clear(self):
        for i in range(self.col): self.table[i].data = list()
        self.row = 0

    def Cov(self, i, j=None):
        if j is None: j = i
        import numpy as np
        xi = np.array(self.table[i].data) - self.table[i].mean()
        xj = np.array(self.table[j].data) - self.table[j].mean()
        return 1. / (self.row - 1) * np.dot(xi.transpose(), xj)

    def R(self, i, j=None):
        from math import sqrt
        if j is None: j = i
        return self.Cov(i, j) / sqrt(self.Cov(i)) / sqrt(self.Cov(j))

    def readCSVFile(self, _fileName, headerstyle=0, _tableName=None, _delimeter=',', typename="double"):
        import re

        def searchNameAndUnitDeli(string, deli=","):
            matchName = re.search(r'^.*?' + deli, string)
            if matchName:
                name = re.search(r'^.*?' + deli, string).group(0)[:-1]
            else:
                name = string.strip()

            matchUnit = re.search(deli + r'.*$', string)
            if matchUnit:
                unit = re.search(deli + r'.*$', string).group(0)[1:]
            else:
                unit = "/"

            return name, unit

        def searchNameAndUnitSqure(string):
            matchName = re.search(r'^.*?\(', string)
            if matchName:
                name = re.search(r'^.*?\(', string).group(0)[:-1]
            else:
                name = string.strip()

            matchUnit = re.search(r'\(.*\)$', string)
            if matchUnit:
                unit = re.search(r'\(.*\)$', string).group(0)[1:-1]
            else:
                unit = "/"

            return name, unit

        self.__init__()
        import csv
        with open(_fileName, 'r') as csvfile:
            content = [each for each in csv.reader(csvfile, delimiter=_delimeter)]

        if _tableName is None:
            row = -1
            for each in content:
                row += 1
                if each[0].strip() != "":
                    rowlabel = row
                    print("Find table in file {} at line {}".format(_fileName, row))
                    break
            if row >= len(content):
                raise Exception("Can not find table in file {}".format(_fileName))

            content[row] = [i for i in content[row] if len(str(i).strip()) != 0]
            self.col = len(content[row])

            try:
                while content[row][0].strip() != "": self.row += 1;row += 1
            except:
                pass
            self.row -= 2
            print("This table has {} rows and {} columns".format(self.row, self.col))

            for i in range(self.col):
                row = rowlabel
                if headerstyle == 0:
                    coltoapp = PhsicalVarList(list(), content[row][i].strip(), content[row + 1][i].strip());
                elif headerstyle == 1:
                    colname, colunit = searchNameAndUnitDeli(content[row][i].strip())
                    coltoapp = PhsicalVarList(list(), colname, colunit);
                elif headerstyle == 2:
                    colname, colunit = searchNameAndUnitSqure(content[row][i].strip())
                    coltoapp = PhsicalVarList(list(), colname, colunit)
                else:
                    raise Exception("Table name and unit delimiter error!!")

                row += 2
                for j in range(self.row):
                    if typename == "double":
                        coltoapp.data.append(eval(content[row + j][i].strip()))
                    elif typename == "string":
                        coltoapp.data.append(content[row + j][i].strip(_delimeter))
                self.table.append(coltoapp)

        else:
            row = 0
            for each in content:
                row += 1
                if _tableName in each:
                    rowlabel = row
                    print("Find table {} in the file {} at line {}".format(_tableName, _fileName, rowlabel))
                    break
            if row >= len(content):
                raise Exception("Can not find the table {} in file {}".format(_tableName, _fileName))

            content[row] = [i for i in content[row] if len(str(i).strip()) != 0]
            self.col = len(content[row])

            try:
                while (content[row][0].strip() != ""): self.row += 1;row += 1
            except:
                pass
            self.row -= 2
            print("This table has {} rows and {} columns".format(self.row, self.col))

            for i in range(self.col):
                row = rowlabel
                coltoapp = PhsicalVarList(list(), content[row][i].strip(), content[row + 1][i].strip());
                row += 2
                for j in range(self.row):
                    coltoapp.data.append(eval(content[row + j][i].strip()))
                self.table.append(coltoapp)
            self.setUnitToSI()

    def readExcelFile(self, filename, sheetname, datatype="string", orient=0):
        from pandas import read_excel
        if orient == 1:
            read_excel("")

    def readSQLiteTable(self, databasename, tablename):
        import sqlite3
        conn = sqlite3.connect(databasename)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(%s)" % tablename)
        header = cursor.fetchall()
        if len(header) == 0:
            raise Exception("No table named {} is in the database {}".format(tablename, databasename))
        print("Table information")
        for row in header:
            print(row)

        name = [each[1] for each in header]  # 读取名字

        cursor.execute("select * from %s" % tablename)
        content = cursor.fetchall()

        print("This table has {} columns and {} rows".format(len(name), len(content)))
        self.clear()
        for i in range(len(name)):
            self.table.append(PhsicalVarList(list()))
        self.col = len(name)
        import re
        def searchNameAndUnitSqure(string):
            matchName = re.search(r'^.*?\(', string)
            if matchName:
                name = re.search(r'^.*?\(', string).group(0)[:-1]
            else:
                name = string.strip()

            matchUnit = re.search(r'\(.*\)$', string)
            if matchUnit:
                unit = re.search(r'\(.*\)$', string).group(0)[1:-1]
            else:
                unit = "/"

            return name, unit

        ColUnit = []
        ColName = []

        for each in name:
            tempname, tempunit = searchNameAndUnitSqure(each)
            ColName.append(tempname)
            ColUnit.append(tempunit)
        # ColUnit = [re.search(r'\(.*\)$', each).group(0)[1:-1] for each in name]
        # ColName = [re.search(r'^.*\(', each).group(0)[:-1] for each in name]
        self.setTableHeader(ColName)
        self.setTableUnit(ColUnit)

        def validvalue(li):
            for each in li:
                if str(each).strip() == "" or each is None:
                    return False
            return True

        # for each in content:
        #     if validvalue(each):
        #         self.append(each)

        for each in content:
            self.append(each)
        conn.close()

    def selectValidData(self, colist=[], isNumber=False):
        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                pass
            return False

        def validvalue(li, isNumber=False):
            if isNumber == False:
                for each in li:
                    if str(each).strip() == "" or (each is None):
                        return False
                return True
            else:
                for each in li:
                    if str(each).strip() == "" or (each is None) or (not is_number(str(each).strip())):
                        return False
                return True

        if len(colist) == 0:
            colist = range(self.col)

        result = ArrayTable(len(colist))
        name, unit = self.getHeader()
        result.setTableHeader([name[j] for j in colist])
        result.setTableUnit([unit[j] for j in colist])

        for i in range(self.row):
            record = self.getOneRecord(i)
            recordpart = [record[j] for j in colist]
            if validvalue(recordpart, isNumber):
                result.append(recordpart)

        return result

    def plot(self, _coly=1, _colx=0):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))

        if type(_coly) == int:
            plt.plot(self.table[_colx].data, self.table[_coly].data)
            plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
            plt.ylabel(self.table[_coly].ColName + "(" + self.table[_coly].ColUnit + ")")
        elif type(_coly) == list:
            for i in _coly:
                plt.plot(self.table[_colx].data, self.table[i].data, label=self.table[i].ColName)
                plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
                plt.legend()
                plt.ylabel(self.table[_coly[0]].ColName + "(" + self.table[_coly[0]].ColUnit + ")")
        # plt.grid()
        # plt.axhline(y=0,color='r',linestyle="-.")
        # plt.axvline(x=0,color='g',linestyle=":")
        plt.tight_layout()
        plt.show()
        return plt

    def scatter(self, _coly=1, _colx=0, hight=None):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))

        if type(_coly) == int:
            plt.scatter(self.table[_colx].data, self.table[_coly].data, s=3)
            if hight is not None:
                plt.scatter(self.table[_colx].data[hight], self.table[_coly].data[hight], s=6, marker="*", c="red")
            plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
            plt.ylabel(self.table[_coly].ColName + "(" + self.table[_coly].ColUnit + ")")
        elif type(_coly) == list:
            for i in _coly:
                plt.scatter(self.table[_colx].data, self.table[i].data, label=self.table[i].ColName)
                plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
                plt.legend()
                plt.ylabel(self.table[_coly[0]].ColName + "(" + self.table[_coly[0]].ColUnit + ")")
        # plt.grid()
        # plt.axhline(y=0,color='r',linestyle="-.")
        # plt.axvline(x=0,color='g',linestyle=":")
        plt.tight_layout()
        plt.show()
        return plt

    def colorScatter(self, _colx=0, _coly=1, _colz=2, hightlight=None):
        import matplotlib.pyplot as plt
        # from matplotlib import cm
        plt.figure(1, figsize=(10, 5))

        plt.scatter(self.table[_colx].data, self.table[_coly].data, s=self.table[_colz].scaleToRange([3, 100]),
                    c='blue')
        if hightlight is not None:
            plt.scatter(self.table[_colx].data[hightlight], self.table[_coly].data[hightlight],
                        s=max(self.table[_colz].scaleToRange([3, 100])),
                        c='red', marker="*")
        plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
        plt.ylabel(self.table[_coly].ColName + "(" + self.table[_coly].ColUnit + ")")
        plt.tight_layout()
        # plt.colorbar()
        plt.show()

    def scatter3D(self, colx=0, coly=1, colz=2):
        import matplotlib.pyplot as plt
        ax = plt.subplot(111, projection='3d')
        for i in range(self.row):
            ax.scatter(self.table[colx].data[i], self.table[coly].data[i], self.table[colz].data[i], c='b')
        plt.tight_layout()
        plt.show()

    def interpolate2D(self, colx=0, coly=1, colz=2):
        import numpy as np
        from scipy import interpolate
        x = np.linspace(self.findMinValue(colx) - 0.2 * abs(self.findMinValue(colx)),
                        self.findMaxValue(colx) + 0.2 * abs(self.findMaxValue(colx)), 100)
        y = np.linspace(self.findMinValue(coly) - 0.2 * abs(self.findMinValue(coly)),
                        self.findMaxValue(coly) + 0.2 * abs(self.findMaxValue(coly)), 100)
        # func=interpolate.interp2d(self.table[colx].data,self.table[coly].data,self.table[colz].data,kind="cubic")
        func = interpolate.Rbf(self.table[colx].data, self.table[coly].data, self.table[colz].data, kind="multiquadric")
        xx, yy = np.meshgrid(x, y)
        zz = func(x, y)
        import matplotlib.pyplot as plt
        # from matplotlib import cm
        fig = plt.figure(figsize=(8, 8))
        plt.contourf(xx, yy, zz, 10, cmap="Blues")
        plt.tight_layout()
        plt.show()

    def plotfun(self, fun, start=0, end=10, step=None):
        from numpy import arange
        if step is None:
            step = (end - start) / 1e3
        result = ArrayTable(2, 0)
        for t in arange(start, end, step):
            result.append([t, fun(t)])

        return result

    def compareResultPlot(self, ArrayList):
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))

        for i in range(1, self.col):
            plt.plot(self.table[0].data, self.table[i].data, label=self.table[i].ColName)
        plt.xlabel(self.table[0].ColName + "(" + self.table[0].ColUnit + ")")
        plt.ylabel(self.table[1].ColName + "(" + self.table[1].ColUnit + ")")

        for each in ArrayList:
            for i in range(1, each.col):
                plt.plot(each.table[0].data, each.table[i].data, label=each.table[i].ColName)

        plt.legend()

        plt.tight_layout()
        plt.show()
        return plt

    def animation(self, _coly=1, _colx=0):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        fig = plt.figure(0, figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)
        # fig, ax = plt.subplots()
        ln, = ax.plot([], [], 'b-', animated=False)
        plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
        plt.ylabel(self.table[_coly].ColName + "(" + self.table[_coly].ColUnit + ")")
        plt.tight_layout()

        def init():
            ax.set_xlim(0.8 * min(self.table[_colx].data), 1.05 * max(self.table[_colx].data))
            ax.set_ylim(0.8 * min(self.table[_coly].data), 1.05 * max(self.table[_coly].data))
            return ln,
            # ax.set_xlim(0, 2 * np.pi)
            # ax.set_ylim(-1, 1)
            # return ln,

        def update(i):
            # ln.set_data(range(i),range(i))
            ln.set_data(self.table[_colx].data[:i], self.table[_coly].data[:i])
            # print(i)
            return ln,

        ani = FuncAnimation(fig, update, frames=np.arange(0, self.row), interval=0.00001,
                            init_func=init, blit=True, repeat=False)
        # ani.save("temp12.gif", writer="pillow",fps=30)
        plt.show()

    def pie(self, _row):
        head, unit = self.getHeader()
        data = self.getOneRecord(_row)
        import matplotlib.pyplot as plt
        plt.figure(1, figsize=(10, 5))
        plt.pie(x=data, labels=head, autopct="%.3f%%")
        plt.show()

    def Log_plot(self, _coly=1, _colx=0):
        import matplotlib.pyplot as plt
        import numpy as np
        plt.figure(1, figsize=(10, 5))
        if type(_coly) == int:
            plt.plot(np.log(self.table[_colx].data), np.log(self.table[_coly].data))
            plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
            plt.ylabel(self.table[_coly].ColName + "(" + self.table[_coly].ColUnit + ")")
        elif type(_coly) == list:
            for i in _coly:
                plt.plot(np.log(self.table[_colx].data), np.log(self.table[i].data), label=self.table[i].ColName)
                plt.xlabel(self.table[_colx].ColName + "(" + self.table[_colx].ColUnit + ")")
                plt.legend()
                plt.ylabel(self.table[_coly[0]].ColName + "(" + self.table[_coly[0]].ColUnit + ")")
        plt.grid()
        plt.show()

    def findMaxValue(self, _col):
        return max(self.table[_col].data)

    def findMinValue(self, _col):
        return min(self.table[_col].data)

    def findMaxValueIndex(self, _col):
        return self.table[_col].data.index(self.findMaxValue(_col))

    def findMinValueIndex(self, _col):
        return self.table[_col].data.index(self.findMinValue(_col))

    def findMaximumDataIndex(self, col, order):
        result = list()
        for j in range(order, self.row - order):
            flag = True
            for i in range(j - order, j + order + 1):
                if self.table[col].data[j] < self.table[col].data[i]:
                    flag = False
                    break
            if flag:
                result.append(j)
        return result

    def getOneRecord(self, _row):
        if _row >= self.row:
            raise Exception(
                "Table only has {} rows but you are trying to get row {} 's data".format(self.row, _row + 1))
        else:
            result = list()
            for i in range(self.col):
                result.append(self.table[i].data[_row])
            return result

    def slice(self, left, right, _col=0):
        if left < min(self.table[_col].data):
            left = min(self.table[_col].data)

        if right > max(self.table[_col].data):
            right = max(self.table[_col].data)

        self.doQuickSort(_col)
        result = ArrayTable(self.col, 0)
        head, unit = self.getHeader()
        result.setTableHeader(head)
        result.setTableUnit(unit)
        start = 0
        while self.table[_col].data[start] < left:
            start += 1

        while self.table[_col].data[start] < right:
            result.append(self.getOneRecord(start))
            start += 1

        return result

    def convertDataToArray(self):
        result = list()
        for i in range(self.row):
            result.append(self.getOneRecord(i))
        return result

    def getHeader(self):
        Name = list()
        Unit = list()
        for i in range(self.col):
            Name.append(self.table[i].ColName)
            Unit.append(self.table[i].ColUnit)
        return Name, Unit

    def writeToCSVFile(self, _fileName):
        import csv
        name, unit = self.getHeader()
        data = self.convertDataToArray()
        with open(_fileName, 'w', newline='') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(name)
            f_csv.writerow(unit)
            f_csv.writerows(data)

    def writeToSQLite(self, filename="temp.db", tablename="temptable"):
        import sqlite3
        conn = sqlite3.connect(filename)
        cursor = conn.cursor()
        cursor.execute("DROP TABLE IF EXISTS `%s`;" % tablename)
        s = "CREATE TABLE `%s`(" % tablename
        for i in range(self.col):
            s += "`{}({})` FLOAT,".format(self.table[i].ColName, self.table[i].ColUnit)
        s = s[:-1] + ");"
        cursor.execute(s)
        for i in range(self.row):
            com = "INSERT INTO %s VALUES (" % tablename
            for j in range(self.col):
                com += "{},".format(self.table[j].data[i])
            com = com[:-1] + ");"
            cursor.execute(com)
        conn.commit()
        conn.close()

    def appendToSQL(self, filename="temp.db", tablename="temptable"):
        import sqlite3
        conn = sqlite3.connect(filename)
        cursor = conn.cursor()
        # cursor.execute("DROP TABLE IF EXISTS `%s`;" % tablename)
        # s = "CREATE TABLE `%s`(" % tablename
        # for i in range(self.col):
        #     s += "`{}({})` FLOAT,".format(self.table[i].ColName, self.table[i].ColUnit)
        # s = s[:-1] + ");"
        # cursor.execute(s)
        for i in range(self.row):
            com = "INSERT INTO %s VALUES (" % tablename
            for j in range(self.col):
                com += "{},".format(self.table[j].data[i])
            com = com[:-1] + ");"
            cursor.execute(com)
        conn.commit()
        conn.close()

    def insertIntoSQLDB(self, databasename, tablename, primarykey=None):
        import sqlite3
        conn = sqlite3.connect(databasename)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(%s)" % tablename)
        header = cursor.fetchall()
        if len(header) == 0:
            raise Exception("No table named {} is in the database {}".format(tablename, databasename))
        print("Table information")
        for row in header:
            print(row)

        databaseHeader = [each[1] for each in header]  # 读取名字

        import re
        def searchNameAndUnitSqure(string):
            matchName = re.search(r'^.*?\(', string)
            if matchName:
                name = re.search(r'^.*?\(', string).group(0)[:-1]
            else:
                name = string.strip()

            matchUnit = re.search(r'\(.*\)$', string)
            if matchUnit:
                unit = re.search(r'\(.*\)$', string).group(0)[1:-1]
            else:
                unit = "/"

            return name, unit

        DBColUnit = []
        DBColName = []
        for each in databaseHeader:
            tempName, tempUnit = searchNameAndUnitSqure(each)
            DBColName.append(tempName)
            DBColUnit.append(tempUnit)

        print(DBColName)
        Name, Unit = self.getHeader()
        print(Name)
        existheaderIndex = []
        for i in range(len(databaseHeader)):
            if DBColName[i] in Name:
                existheaderIndex.append(i)
        existheaderIndexInTable = [Name.index(i) for i in [DBColName[j] for j in existheaderIndex]]

        existheader = "(`" + "`,`".join([databaseHeader[i] for i in existheaderIndex]) + "`)"

        print("Following columns are included the database:", existheader)

        if primarykey is None:
            for i in range(self.row):
                command = "INSERT INTO " + tablename + existheader + "values" + "("
                for j in existheaderIndexInTable:
                    command += '"' + str(self.table[j].data[i]) + '",'
                command = command[:-1] + ")"
                cursor.execute(command)

        conn.commit()
        conn.close()

    def setTableHeader(self, header):
        if len(header) != self.col:
            raise Exception("Header's size does not match the table size")
        else:
            for i in range(self.col):
                self.table[i].ColName = header[i]

    def setTableUnit(self, unitHeader):
        if len(unitHeader) != self.col:
            raise Exception("Unit header's size does not match the table size")
        else:
            for i in range(self.col):
                self.table[i].ColUnit = unitHeader[i]

    def integrate(self, _coly, end_row=None, start_row=0, _colx=0) -> float:
        if end_row is None: end_row = self.row - 1
        if end_row > self.row - 1:
            raise Exception(
                "can not integrate to row {} because this table only has {} rows".format(end_row, self.row - 1))
        result = 0.
        for i in range(start_row, end_row):
            result += (self.table[_colx].data[i + 1] - self.table[_colx].data[i]) * (
                    self.table[_coly].data[i + 1] + self.table[_coly].data[i]) / 2.
        return result

    def linearInterpolate(self, x, col, type=0):
        left = 0
        right = self.row - 1
        while left != right - 1:
            mid = (left + right) // 2
            if self.table[0].data[mid] <= x:
                left = mid
            elif self.table[0].data[mid] >= x:
                right = mid

        x1 = self.table[0].data[left]
        x2 = self.table[0].data[right]
        y1 = self.table[col].data[left]
        y2 = self.table[col].data[right]
        return y1 * (x - x2) / (x1 - x2) + y2 * (x - x1) / (x2 - x1)

    def selectColumns(self, colist):
        result = ArrayTable(len(colist), 0)
        head, unit = self.getHeader()
        headS = [head[i] for i in colist]
        units = [unit[i] for i in colist]
        result.setTableHeader(headS)
        result.setTableUnit(units)
        for i in range(len(colist)):
            result.table[i].data = self.table[colist[i]].data
        result.col = len(colist)
        result.row = self.row
        return result

    def polyFit(self, coly=1, colx=0, _deg=1):
        from numpy import polyfit, polyval
        para = polyfit(self.table[colx].data, self.table[coly].data, deg=_deg)
        for i in range(len(para)):
            print("({}*x^{})+".format(para[i], _deg - i), end="")
        print("\b")

        result = self.selectColumns([colx, coly])
        colToApp = PhsicalVarList([], "Fitted data", result.table[1].ColUnit)
        colToApp.data = polyval(para, result.table[0].data)
        result.appendColumn(colToApp)
        return result, para

    def createSimilarEmptyTable(self, colist):
        result = ArrayTable(len(colist), 0)
        head, unit = self.getHeader()
        headS = [head[i] for i in colist]
        units = [unit[i] for i in colist]
        result.setTableHeader(headS)
        result.setTableUnit(units)
        return result

    def setUnitToSI(self):
        for i in range(self.col):
            if self.table[i].ColUnit.strip() == "mm":
                self.table[i].ColUnit = "m"
                self.table[i] * 1.e-3
            if self.table[i].ColUnit.strip() == "bar":
                self.table[i].ColUnit = "Pa"
                self.table[i] * 1.e5
            if self.table[i].ColUnit.strip() == "MPa":
                self.table[i].ColUnit = "bar"
                self.table[i] * 1.e1

    def diff(self, _coly=1, _colx=0):
        temp = PhsicalVarList([0.] * self.row)
        temp.ColName = "D(" + self.table[_coly].ColName + ")D(" + self.table[_colx].ColName + ")"
        temp.ColUnit = self.table[_coly].ColUnit + "/" + self.table[_colx].ColUnit;
        temp.data[0] = (self.table[_coly].data[1] - self.table[_coly].data[0]) / (
                self.table[_colx].data[1] - self.table[_colx].data[0])
        for i in range(1, self.row - 1):
            temp.data[i] = (self.table[_coly].data[i + 1] - self.table[_coly].data[i - 1]) / (
                    self.table[_colx].data[i + 1] - self.table[_colx].data[i - 1])
        temp.data[-1] = (self.table[_coly].data[- 1] - self.table[_coly].data[- 2]) / (
                self.table[_colx].data[- 1] - self.table[_colx].data[- 2])
        self.table.append(temp)
        self.col += 1

    def delRow(self, _row):
        if _row >= self.row:
            raise Exception("can't delete the row because it does not exist!!")
        result = self.getOneRecord(_row)
        for i in range(self.col):
            del self.table[i].data[_row]
        self.row -= 1
        return result

    def popRow(self):
        # self.row -= 1
        # return self.delRow(self.row - 1)

        result = self.getOneRecord(self.row - 1)
        for i in range(self.col):
            self.table[i].data.pop()
        self.row -= 1
        return result

    def openWithProgram(self):
        self.writeToCSVFile("temp.csv")
        from os import system
        system("start temp.csv")
        system("pause")
        system("del /s /q temp.csv")

    def translate(self, jsonfile="NameDic.json"):
        def getKey(dic, value):
            if value not in dic.values():
                return None
            else:
                for i in dic.keys():
                    if dic[i] == value: return i

        def getKey2(dic, value):
            for each in dic.values():
                if value in each:
                    return getKey(dic, each)
            else:
                return None

        import json
        with open(jsonfile, 'r+', encoding='utf-8') as f:
            content = json.load(f)

        for i in range(self.col):
            key = getKey2(content, self.table[i].ColName)
            if key is not None:
                self.table[i].ColName = key

    def exchangeRow(self, i, j):
        temp = self.getOneRecord(i)
        for k in range(len(temp)):
            self.table[k].data[i] = self.table[k].data[j]

        for k in range(len(temp)):
            self.table[k].data[j] = temp[k]

    def doQuickSort(self, col=0):

        self.__quickSort(0, self.row - 1, col)

    def __quickSort(self, left, right, col):
        if left >= right: return
        i = left
        j = right
        base = self.getOneRecord(left)

        while i < j:
            while self.table[col].data[j] >= base[col] and i < j:
                j -= 1
            while self.table[col].data[i] <= base[col] and i < j:
                i += 1
            if i < j:
                self.exchangeRow(i, j)

        for k in range(len(base)):
            self.table[k].data[left] = self.table[k].data[i]

        for k in range(len(base)):
            self.table[k].data[i] = base[k]

        self.__quickSort(left, i - 1, col)
        self.__quickSort(i + 1, right, col)
