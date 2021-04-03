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
        self.table.append(_data)
        self.col += 1

    def clear(self):
        for i in range(self.col): self.table[i].data = list()
        self.row = 0

    def readCSVFile(self, _fileName, _tableName=None, _delimeter=','):
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

    def readSQLiteTable(self, databasename, tablename):
        import sqlite3
        conn = sqlite3.connect(databasename)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(%s)" % tablename)
        header = cursor.fetchall()
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
        self.setTableHeader(name)

        for each in content:
            self.append(each)
        conn.close()

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
        plt.show()

    def animation(self, _coly=1, _colx=0):
        pass

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

    def animation(self):
        pass

    def findMaxValue(self, _col):
        return max(self.table[_col].data)

    def findMinValue(self, _col):
        return min(self.table[_col].data)

    def findMaxValueIndex(self, _col):
        return self.table[_col].data.index(self.findMaxValue(_col))

    def findMinValueIndex(self, _col):
        return self.table[_col].data.index(self.findMinValue(_col))

    def getOneRecord(self, _row):
        if _row >= self.row:
            raise Exception(
                "Table only has {} rows but you are trying to get row {} 's data".format(self.row, _row + 1))
        else:
            result = list()
            for i in range(self.col):
                result.append(self.table[i].data[_row])
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

    def integral(self, _coly, _colx=0) -> float:
        result = 0.
        for i in range(self.row - 1):
            result += (self.table[_colx].data[i + 1] - self.table[_colx].data[i]) * (
                    self.table[_coly].data[i + 1] - self.table[_coly].data[i]) / 2.
        return result

    def setUnitToSI(self):
        for i in range(self.col):
            if self.table[i].ColUnit.strip() == "mm":
                self.table[i].ColUnit = "m"
                self.table[i] * 1.e-3

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

