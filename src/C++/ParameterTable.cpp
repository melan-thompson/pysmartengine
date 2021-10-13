#include "ParameterTable.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include "Functions.h"

#define PI 3.141592654
using namespace function;

bool ParameterTable::append(std::string paraName, std::string paraUnit, double value)
{
    num_of_row++;
    table.push_back(para(paraName, paraUnit, value));
    return true;
}

////////////////////////Definitions of ParameterTable functions/////////////////////////////////////////////

bool ParameterTable::readInFile(std::string filename, std::string start, int row)
{
    using namespace std;
    string original;

    ifstream fin;
    fin.open(filename);
    if (!fin.is_open())
    {
        cerr << "Trying to open file " << filename << ", but fail to open that file" << endl;
        system("pause");
        exit(EXIT_FAILURE);
    }

    if (start != "")//读取表格名称
    {
        while (fin.good())
        {
            getline(fin, original, '\n');
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
        }
    }

    if (row == 0)
    {
        while (fin.good())
        {
            getline(fin, original, '\n');
            if (!(Trim(original) == ""))
            {
                row++;
            }
            else
            {
                break;
            }
        }
    }

    num_of_row = row;
    table.resize(num_of_row);

    fin.close();
    fin.open(filename);
    if (start != "")
    {
        while (fin.good())
        {
            getline(fin, original, '\n');
            if (Trim(original) == start)
            {
                break;
            }
        }
    }

    for (int i = 0; i < row; i++)
    {

        getline(fin, original, ',');
        table[i].parameter_name = Trim(original);
        getline(fin, original, ',');
        table[i].parameter_unit = Trim(original);
        getline(fin, original, '\n');
        if (Trim(original) == "")
        {
            cerr << table[i].parameter_name << " parameter value is not given, will return -inf!!!" << endl;
            system("pause");
            table[i].parameter_value = -INFINITY;
        }
        else if (!isNumber(Trim(original)))
        {
            cerr << table[i].parameter_name << " parameter value is not a number, will return -inf!!!" << endl;
            system("pause");
            table[i].parameter_value = -INFINITY;
        }
        else
        {
            stringstream ss;
            ss << Trim(original);
            ss >> table[i].parameter_value;
        }
    }
    fin.close();
    setUnit2SI();
    return true;
}

ParameterTable::ParameterTable(std::string filename, std::string start, int row)
{
    readInFile(filename, start, row);
}

ParameterTable::ParameterTable(int row, std::string name)
{
    num_of_row = row;
    table.resize(row);
    table_name = name;
}

template <class T>
bool ParameterTable::setParasName(T& str)
{
    for (int i = 0; i < getArrayLen(str) && i < num_of_row; i++)
    {
        table[i].parameter_name = str[i];
    }
    return true;
}


template <class T>
bool ParameterTable::setParasUnit(T& str)
{
    for (int i = 0; i < getArrayLen(str) && i < num_of_row; i++)
    {
        table[i].parameter_unit = str[i];
    }
    return true;
}

template <class T>
bool ParameterTable::setParasValue(T& arr)
{
    for (int i = 0; i < getArrayLen(arr) && i < num_of_row; i++)
    {
        table[i].parameter_value = arr[i];
    }
    return true;
}

bool ParameterTable::writeToFile(std::string filename, int mode, std::string title)
{
    using namespace std;
    ofstream fout;
    if (mode == 0)//重写文件
    {
        fout.open(filename, ios_base::out | ios_base::trunc);
    }
    else if (mode == 1)//文件后追加
    {
        fout.open(filename, ios_base::out | ios_base::app);
    }

    if (!fout.is_open())
    {
        cerr << "Fail to open the output file: " << filename << endl;
        system("pause");
        exit(EXIT_FAILURE);
    }

    if (title != "")
    {
        fout << title << ",," << "\n";
    }
    //fout << table_name << ",," << "\n";
    for (int i = 0; i < num_of_row; i++)
    {
        fout << table[i].parameter_name << "," << table[i].parameter_unit << "," << table[i].parameter_value << "\n";
    }

    fout << "\n";
    return true;
}


void ParameterTable::setUnit2SI()
{
    for (int i = 0; i < num_of_row; i++)
    {
        if (table[i].parameter_unit == "℃")
        {
            table[i].parameter_unit = 'K';
            table[i].parameter_value += 273.15;
        }
        else if (table[i].parameter_unit == "MPa")
        {
            table[i].parameter_unit = "Pa";
            table[i].parameter_value *= pow(10, 6);
        }
        else if (table[i].parameter_unit == "bar")
        {
            table[i].parameter_unit = "Pa";
            table[i].parameter_value *= pow(10, 5);
        }
        else if (table[i].parameter_unit == "mm")
        {
            table[i].parameter_unit = "m";
            table[i].parameter_value *= 1e-3;
        }
        else if (table[i].parameter_unit == "kW")
        {
            table[i].parameter_unit = "W";
            table[i].parameter_value *= 1e3;
        }
        else if (table[i].parameter_unit == "mm^3")
        {
            table[i].parameter_unit = "m^3";
            table[i].parameter_value *= 1e-9;
        }
        else if (table[i].parameter_unit == "mg")
        {
            table[i].parameter_unit = "kg";
            table[i].parameter_value *= 1e-6;
        }
    }
}

void ParameterTable::show()
{
    using namespace std;
    cout << "           -------------------------------------------------------------               " << endl;
    cout << setw(46) << table_name << endl;
    cout << "           -------------------------------------------------------------               " << endl;
    for (int i = 0; i < num_of_row; i++)
    {
        cout << setw(30) << table[i].parameter_name << setw(20) << table[i].parameter_unit << setw(20) << setiosflags(ios::fixed) << setprecision(5) << table[i].parameter_value << endl;
    }
    cout << "           -------------------------------------------------------------               " << endl;
}

double ParameterTable::searchValue(std::string index)
{
    for (int i = 0; i < num_of_row; i++)
    {
        if (table[i].parameter_name == index) {
            return table[i].parameter_value;
        }
    }
    std::cout << index << " parameter does not exist!!! Please check the input file !!!" << std::endl;
    system("pause");
    exit(EXIT_FAILURE);
    return -INFINITY;
}

std::string ParameterTable::searchUnit(std::string index)
{
    for (int i = 0; i < num_of_row; i++)
    {
        if (table[i].parameter_name == index) {
            return table[i].parameter_unit;
        }
    }
    std::cout << index << " parameter does not exist!!! Please check the input file !!!" << std::endl;
    system("pause");
    exit(EXIT_FAILURE);
    return "parameter does not exist";
}



double ParameterTable::sum(int start, int end)
{
    double s = 0;
    for (int i = start; i <= end; i++)
    {
        s += table[i].parameter_value;
    }
    return s;
}
////////////////////////////////////////End/////////////////////////////////////////////////////

