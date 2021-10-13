#include "ArrayTable.h"
#define PI 3.141592654
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include "Functions.h"
using namespace function;

template <class Type>
PhysicalVarList<Type>::PhysicalVarList(int row)
{
    ColName = "EmptyName";
    ColUnit = "EmptyUnit";
    if (row != 0)
    {
        for (int i = 1; i <= row; i++)
        {
            data.push_back(0);
        }
    }
}

template <class Type>
bool PhysicalVarList<Type>::add(Type num)
{
    for (int i = 0; i < data.size(); i++)
    {
        data[i] += num;
    }
    return true;
}

template <class Type>
bool PhysicalVarList<Type>::Multiply(double num, std::string new_unit)
{
    for (int i = 0; i < data.size(); i++)
    {
        data[i] *= num;
    }
    ColUnit = new_unit;
    return true;
}

template <class Type>
std::ostream& operator<<(std::ostream& out, PhysicalVarList<Type> L)
{
    out << L.ColName << "\n" << L.ColUnit << "\n";
    for (int i = 0; i < L.data.size(); i++)
    {
        out << L.data[i] << "\n";
    }
    return out;
}

/////////////////////////////////Definitions of ArrayTable functions///////////////////////////
bool ArrayTable::readInFile(std::string filename, std::string start,char delimiter)
{
    using namespace std;
    using namespace function;
    
    string flag;
    string original;

    ifstream fin;
    fin.open(filename,std::ios_base::in);
    if (!fin.is_open())
    {
        cerr << "Trying to open file " << filename << ",but fail to open that file" << endl;
        system("pause");
        exit(EXIT_FAILURE);
    }

    long pos;
    int linenum = 0;
    if (start != "")//读取表格名称
    {
        while (fin.good()&&!fin.eof())
        {
            getline(fin, original, '\n');
            linenum++;
            if (Trim(original) == start)
            {
                break;
            }
        }
        if (fin.peek() == EOF)
        {
            cerr << start << " table does not exist!!" << endl;
            system("pause");
            exit(EXIT_FAILURE);
            return false;
        }
        else
        {
            table_name = Trim(original); 
            std::cout << "Table " << table_name << " found at line:" << linenum << ",byte position:"<< fin.tellg();
        }
    }
    pos = fin.tellg();//记录位置
    
    while (fin.good() && !fin.eof())
    {
        getline(fin, original, '\n');
        if (!(Trim(original) == ""))
        {
            row++;
            static int temp= calStringContainsChar(original, delimiter);
            if (temp!= calStringContainsChar(original, delimiter))
            {
                std::cerr << "This table have some problems, row 1 has " << temp << " columns, but row " << row << " has " << calStringContainsChar(original, delimiter) << " columns!!!\n";
                system("pause");
            }
            col = temp+1;
            
        }
        else
        {
            break;
        }
    }
    row -= 2;

    table.clear();
    table.resize(col);
    std::cout << ". And it has "<<col<<" columns, " << row << " rows to read.\n";
    for (int i = 0; i < col; i++)
    {
        table[i].data.resize(row);
    }
    

    if (fin.tellg() == -1)
    {
        fin.close();
        fin.open(filename, std::ios_base::in);
        fin.seekg(pos, std::ios::beg);
    }

    ////////////读取内容//////////////
    for (int i = 0; i < col-1; i++)
    {
        getline(fin, original, delimiter);
        table[i].ColName = Trim(original);
    }
    getline(fin, original, '\n');
    table[col-1].ColName = Trim(original);

    for (int i = 0; i < col-1; i++)
    {
        getline(fin, original, delimiter);
        table[i].ColUnit = Trim(original);
    }
    getline(fin, original, '\n');
    table[col - 1].ColUnit = Trim(original);


    for (int j = 0; j < row; j++)
    {
        for (int i = 0; i < col-1; i++)
        {
            getline(fin, original, delimiter);
            std::string temp = Trim(original);
            if (temp == "")
            {
                cerr << "Data at row " << j << " column " << i << " is blank, will return -inf!!!" << endl;
                system("pause");
                table[i].data[j] = -INFINITY;
            }
            else if (!isNumber(temp))
            {
                cerr << "Data at row " << j << " column " << i << " is: " << temp << " which is not a number, will return -inf!!!" << endl;
                system("pause");
                table[i].data[j] = -INFINITY;
            }
            else
            {
                table[i].data[j]=str2double(temp);
            }
        }

        getline(fin, original, '\n');
        std::string temp = Trim(original);
        if (temp == "")
        {
            cerr << "Data at row " << j << " column " << col-1 << " is blank, will return -inf!!!" << endl;
            system("pause");
            table[col-1].data[j] = -INFINITY;
        }
        else if (!isNumber(temp))
        {
            cerr << "Data at row " << j << " column " << col-1 << " is: " << temp << " which is not a number, will return -inf!!!" << endl;
            system("pause");
            table[col-1].data[j] = -INFINITY;
        }
        else
        {
            table[col-1].data[j] = str2double(temp);
        }
    }
    fin.close();
    show();
    setUnit2SI();
    return true;
}

ArrayTable::ArrayTable(int _col, int _row, std::string _tablename)
{
    table_name = _tablename;
    table.resize(_col);
    col = _col;
    row = _row;
    for (int i = 0; i < col; i++)
    {
        table[i].data.resize(_row);
    }
}

ArrayTable::ArrayTable(std::string filename, std::string start)
{
    readInFile(filename,start);
}

////////////////////////////以下为表格输出操作///////////////////////////////////
void ArrayTable::show()
{
    using namespace std;
    for (int i = 0; i < col; i++)
    {
        cout << setw(30) << table[i].ColName;
    }
    cout << endl;

    for (int i = 0; i < col; i++)
    {
        cout << setw(30) << table[i].ColUnit;
    }
    cout << endl;

    for (int j = 0; j < row; j++)
    {
        for (int i = 0; i < col; i++)
        {
            cout << setw(30) << table[i].data[j];
        }
        cout << endl;
    }
}

void ArrayTable::show(int column)
{
    using namespace std;
    cout << table[column - 1].ColName << endl;
    cout << table[column - 1].ColUnit << endl;
    for (int i = 0; i < row; i++)
    {
        cout << table[column - 1].data[i] << endl;
    }
}

bool ArrayTable::writeToFile(std::string filename, write_to_file_type mode)
{
    using namespace std;
    using namespace function;
    ofstream fout;
    if (mode == write_to_file_type::truncate2file)//重写文件
    {
        fout.open(filename, ios_base::out | ios_base::trunc /*|ios_base::binary*/);
    }
    else if (mode == write_to_file_type::append2file)//文件后追加
    {
        fout.open(filename, ios_base::out | ios_base::app /*| ios_base::binary*/);
    }

    if (!fout.is_open())
    {
        cerr << "文件" << filename << "打开失败！";
        return false;
        system("pause");
        exit(EXIT_FAILURE);
    }
    //fout << table_name << ",," << "\n";
    for (int i = 0; i < col - 1; i++)
    {

        fout << table[i].ColName + ",";
    }
    fout << table[col - 1].ColName + "\n";

    for (int i = 0; i < col - 1; i++)
    {
        fout << table[i].ColUnit + ",";
    }
    fout << table[col - 1].ColUnit + "\n";

    for (int j = 0; j < row; j++)
    {
        for (int i = 0; i < col - 1; i++)
        {
            fout << table[i].data[j] << ",";
        }
        fout << table[col - 1].data[j] << "\n";

    }
    fout.close();

    return true;
}

bool ArrayTable::plot(int _coly, int _colx)
{
    if (_coly >= col || _colx >= col)
    {
        std::cerr << "Can not plot because data column is out of range," << "table have " << col << " columns," << "but you want to plot at column " << _coly << " to " << _colx << "\n";
        system("pause");
        exit(EXIT_FAILURE);
    }

    using namespace std;
    ofstream fout;
    fout.open(".\\temp.py", ios_base::out | ios_base::trunc);
    string xlist = "["; string ylist = "[";
    for (int i = 0; i < row; i++)
    {
        xlist += to_string(table[_colx].data[i]) + ",";
        ylist += to_string(table[_coly].data[i]) + ",";
    }
    xlist += "]"; ylist += "]";
    string quotation = "\'";
    string str = "#-*- coding:utf-8 -*-\n\n";
    str += "import matplotlib.pyplot as plt\n";
    str += "plt.plot(" + xlist + "," + ylist + ")\n";
    str += "plt.xlabel(" + quotation + table[_colx].ColName + "(" + table[_colx].ColUnit + ")" + quotation + ")\n";
    str += "plt.ylabel(" + quotation + table[_coly].ColName + "(" + table[_coly].ColUnit + ")" + quotation + ")\n";
    str += "plt.grid()\n";
    str += "plt.show()";

    fout << str;
    fout.close();

    {
        fout.open(".\\ANSI2utf8.py", ios_base::out | ios_base::trunc);
        string toUtf8 = "import codecs\n";
        toUtf8 += "import os\n";
        toUtf8 += "f = codecs.open(" + quotation + "temp.py" + quotation + ", 'r', 'ansi')\n";
        toUtf8 += "ff = f.read()\n";
        toUtf8 += "file_object = codecs.open(" + quotation + "temp.py" + quotation + ", 'w', 'utf-8')\n";
        toUtf8 += "file_object.write(ff)\n";
        fout << toUtf8;
        fout.close();
        system("python .\\ANSI2utf8.py");
        system("del .\\ANSI2utf8.py");
    }
    system("python .\\temp.py");
    //system("pause");
    system("del .\\temp.py");
    return true;
}


void ArrayTable::openWithProgram()
{
    writeToFile("TempFile.csv");
    system("start TempFile.csv");
    std::cout << "The table has been opened in your program, please close it first and press enter to continue!\n";
    system("pause");
    system("del /q TempFile.csv");
}
///////////////////////////以上为表格输出操作/////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////表格修改操作//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void ArrayTable::setUnit2SI()
{
    for (int i = 0; i < col; i++)
    {
        if (table[i].ColUnit == "mm")
        {
            table[i].ColUnit = "m";
            for (int j = 0; j < row; j++)
            {
                table[i].data[j] *= 1e-3;
            }
        }
        if (table[i].ColUnit == "bar")
        {
            table[i].ColUnit = "Pa";
            for (int j = 0; j < row; j++)
            {
                table[i].data[j] *= 1.e5;
            }
        }

    }
}

bool ArrayTable::append(std::vector<double> arr)
{
    if (arr.size() != col)
    {
        std::cerr << "The array to be appended does not match the table's size!!! Please check!!" << std::endl;
        system("pause");
        return false;
    }
    row++;
    for (int i = 0; i < col; i++)
    {
        table[i].data.push_back(arr[i]);
    }
    return true;
}

bool ArrayTable::setRowValue(int _row, std::vector<double> _arr)
{
    if (_row > row)
    {
        std::cerr << "the row you are going to set values is out of table range, table have " << row << " rows, but you are going to change row " << _row << "'s value!!\n";
        system("pause");
        return false;
    }
    if (_arr.size() != col)
    {
        std::cerr << "Array size does not match table column, table have " << col << " columns, but the array size is " << _arr.size() << "!!\n";
        system("pause");
        return false;
    }
    for (int i = 0; i < col; i++)
    {
        table[i].data[_row] = _arr[i];
    }
    return true;
}

bool ArrayTable::exchageRowValue(int _row1, int _row2)
{
    if (_row1 > row || _row2 > row)
    {
        std::cerr << "table have " << row << " rows, but you are going the exchage row " << _row1 << " and row " << _row2 << " data!!\n";
        system("pause");
        return false;
    }
    std::vector<double> temp1 = selectRowValue(_row1);
    std::vector<double> temp2 = selectRowValue(_row2);
    setRowValue(_row1, temp2);
    setRowValue(_row2, temp1);
    return true;
}

bool ArrayTable::setParaName(std::vector < std::string> _arr)
{
    for (int i = 0; i < _arr.size() && i < col; i++)
    {
        table[i].ColName = _arr[i];
    }
    return true;
}

bool ArrayTable::setParaUnit(std::vector < std::string> _arr)
{
    for (int i = 0; i < _arr.size() && i < col; i++)
    {
        table[i].ColUnit = _arr[i];
    }
    return true;
}
//////////////////////////以上为表格修改操作//////////////////////////////////////

double ArrayTable::linearInterpolate(double x, int columny, int columnx, int type)
{
    //插值前确保已经排好序
    if (type == 1)
    {
        if (x < table[0].data[0] || x > table[0].data[row - 1]) return 0;
    }

    //循环插值需要做一定处理
    double temp = table[0].data[row - 1] - table[0].data[0];
    if (x < table[0].data[0])
    {
        double temp2 = fmod((table[0].data[0] - x), temp);
        x = table[0].data[row - 1] - temp2;

    }
    else if (x > table[0].data[row - 1])
    {
        double temp2 = fmod((x - table[0].data[row - 1]), temp);
        x = table[0].data[0] + temp2;
    }


    int left = 0;
    int right = row - 1;
    while (left != (right - 1))
    {
        int mid = (left + right) / 2;
        if (table[columnx].data[mid] <= x)
        {
            left = mid;
        }
        else if (table[columnx].data[mid] >= x)
        {
            right = mid;
        }

    }

    double x1 = table[columnx].data[left]; double x2 = table[columnx].data[right];
    double y1 = table[columny].data[left]; double y2 = table[columny].data[right];
    return y1 * (x - x2) / (x1 - x2) + y2 * (x - x1) / (x2 - x1);
}

//double ArrayTable::triInterpolate(double _x, int _coly, int _colx)
//{
//    //这里假设表格已经排好序
//    Eigen::VectorXcd y(row);
//    for (int i = 0; i < row; i++)
//    {
//        y(i) = table[_coly].data[i];
//    }
//    Eigen::VectorXcd ffty = fft(y);
//    double xtemp = (_x - table[_colx].data[0]) / double(table[_colx].data[row - 1] - table[_colx].data[0]);
//    int n = row + 1; double result = 0;
//    if (n%2==0)
//    {
//        result += ffty(0).real() + ffty(n / 2).real() * cos(xtemp * n * PI);
//        for (int i = 1; i <= n/2-1; i++)
//        {
//            result += 2. * (ffty(i).real() * cos(xtemp * 2 * PI * i) - ffty(i).imag() * sin(xtemp * 2 * PI * i));
//        }
//    }
//    else
//    {
//        result += ffty(0).real();
//        for (int i = 1; i <= (n-1)/2; i++)
//        {
//            result += 2. * (ffty(i).real() * cos(xtemp * 2 * PI * i) - ffty(i).imag() * sin(xtemp * 2 * PI * i));
//        }
//    }
//    return result/n;
//}

double ArrayTable::Integrate(int _coly, int _start, int _end, int _colx)
{
    if (_start < 0)
    {
        _start = 0;
    }
    if (_end >= row)
    {
        _end = row - 1;
    }
    double Summ = 0;
    for (int i = _start; i < _end - 1; i++)
    {
        Summ += (table[_coly].data[i] + table[_coly].data[i + 1]) / 2 * (table[_colx].data[i + 1] - table[_colx].data[i]);
    }
    return Summ;
}

void ArrayTable::diff(int _coly, int _colx)
{
    PhysicalVarList<double> temp(row);
    temp.ColName = "D(" + table[_coly].ColName + ")D(" + table[_colx].ColName + ")";
    temp.ColUnit = table[_coly].ColUnit + "/" + table[_colx].ColUnit;
    temp.data[0] = (table[_coly].data[1] - table[_coly].data[0]) / (table[_colx].data[1] - table[_colx].data[0]);
    for (int i = 1; i < row - 1; i++)
    {
        temp.data[i] = (table[_coly].data[i + 1] - table[_coly].data[i - 1]) / (table[_colx].data[i + 1] - table[_colx].data[i - 1]);
    }
    temp.data[row - 1] = (table[_coly].data[row - 1] - table[_coly].data[row - 2]) / (table[_colx].data[row - 1] - table[_colx].data[row - 2]);
    table.push_back(temp);
    col += 1;
}



////////////////////////////////以下为表格查询操作////////////////////////////////////////
std::vector<int> ArrayTable::findMaximumDataIndex(int _col, int _order)
{
    std::vector<int> result;
    for (int j = _order; j < row - _order; j++)
    {
        bool flag = true;
        for (int i = j - _order; i <= j + _order; i++)
        {
            if (table[_col].data[j] < table[_col].data[i])
            {
                flag = false; break;
            }
        }
        if (flag)
        {
            result.push_back(j);
        }
    }
    return result;
}

std::vector<double> ArrayTable::selectRowValue(int _row)//读取行row的一条记录
{
    std::vector<double> result(col);
    if (_row >= row)
    {
        std::cerr << "You can't get row " << _row << " record, because table only have " << row << " rows!\n";
        std::cerr << "Will return vector with random result!!!\n";
        system("pause");
        return result;
    }
    for (int i = 0; i < col; i++)
    {
        result[i] = table[i].data[_row];
    }
    return result;
}

void ArrayTable::pop()
{
    row -= 1;
    for (int i = 0; i < col; i++)
    {
        table[col].data.pop_back();
    }
}


void ArrayTable::Clear()
{
    row = 0;
    for (int i = 0; i < col; i++)
    {
        table[i].data.resize(0);
    }
}

int ArrayTable::findMaxDataIndex(int _col)
{
    double temp = 0;
    for (int i = 1; i < row; i++)
    {
        if (table[_col].data[i] > table[_col].data[temp]) temp = i;
    }
    return temp;
}

int ArrayTable::findMinDataIndex(int _col)
{
    double temp = 0;
    for (int i = 1; i < row; i++)
    {
        if (table[_col].data[i] < table[_col].data[temp]) temp = i;
    }
    return temp;
}

int ArrayTable::findMaxDataIndex(int _col) const
{
    double temp = 0;
    for (int i = 1; i < row; i++)
    {
        if (table[_col].data[i] > table[_col].data[temp]) temp = i;
    }
    return temp;
}

int ArrayTable::findMinDataIndex(int _col) const
{
    double temp = 0;
    for (int i = 1; i < row; i++)
    {
        if (table[_col].data[i] < table[_col].data[temp]) temp = i;
    }
    return temp;
}

double ArrayTable::findMaxData(int _col) const
{
    return table[_col].data[findMaxDataIndex(_col)];
}

double ArrayTable::findMinData(int _col) const
{
    return table[_col].data[findMinDataIndex(_col)];
}


double ArrayTable::sumAll(int column) const
{
    double sum = 0;
    for (int i = 0; i < row; i++)
    {
        sum += table[column].data[i];
    }
    return sum;
}

double ArrayTable::quasierror(int _col, double mean) const
{
    double result = 0;
    for (int i = 0; i < row; i++)
    {
        result += pow((table[_col].data[i] - mean), 2);
    }
    return result;
}

void ArrayTable::doQuickSort(int _col)//对表格快速排序
{
    if (_col > col)
    {
        std::cerr << "Can not do QuickSort accordging to column " << _col << ", because table only have " << col << " columns!!\n";
        system("pause");
        return;
    }
    QuickSort(_col, this, 0, row - 1);
}

void ArrayTable::QuickSort(int _col, ArrayTable point[], int low, int high)
{
    if (low < high)
    {
        int i = low;
        int j = high;
        std::vector<double> key = point->selectRowValue(i);
        while (i < j)
        {
            while (i < j && point->table[_col].data[j] >= key[_col])
            {
                j--;
            }
            if (i < j)
            {
                std::vector<double> temp = point->selectRowValue(j);
                point->setRowValue(i, temp);
            }
            while (i < j && point->table[_col].data[i] <= key[_col])
            {
                i++;
            }
            if (i < j)
            {
                std::vector<double> temp = point->selectRowValue(i);
                point->setRowValue(j, temp);
            }
        }
        point->setRowValue(i, key);
        QuickSort(_col, point, low, i - 1);
        QuickSort(_col, point, i + 1, high);
    }
}

//////////////////////////////////End///////////////////////////////////////////////////////
