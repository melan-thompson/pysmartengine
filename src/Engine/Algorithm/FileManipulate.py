def get_pakage_dir(pakagename):
    import sys
    dirname=sys.executable.split("\\")[0:-1]
    dirname.append("Lib")
    dirname.append("site-packages")
    dirname.append(pakagename)
    result="\\".join(dirname)
    import os
    if os.path.isdir(result):
        return result
    else:
        raise Exception("directory does not exist!!")


def TransposeCSVFile(filename):
        def trimEmptylines(content):#去除空行
            rest=[]
            for i in range(len(content)):
                allemptyflag=True
                for j in range(len(content[i])):
                    if content[i][j]!= "":
                        allemptyflag=False
                        break
                if not allemptyflag:
                    rest.append(content[i])
            return rest

        import csv
        filecontent=[]
        with open(filename,'r') as csvFile:
            reader=csv.reader(csvFile)
            for line in reader:
                filecontent.append(line)

        filecontent=trimEmptylines(filecontent)

        result = [[row[i].strip() for row in filecontent] for i in range(len(filecontent[0]))]
        result=trimEmptylines(result)
        print("this table has {} rows {} columns".format(len(result[0]),len(result)))


        with open(filename,'w',newline='') as csvFile:
            writer=csv.writer(csvFile)
            writer.writerows(result)
        print("Successfully transpose file {}".format(filename))