#pragma once
#ifndef PARAMETERTABLE_H
#define PARAMETERTABLE_H
#include <iostream>
#include <vector>



class para
{
public:
    std::string parameter_name;
    std::string parameter_unit;
    double parameter_value;
    para() 
    {
        parameter_name = "EmptyName";
        parameter_unit = "EmptyUnit";
        parameter_value = -INFINITY;
    };
    para(std::string paraName, std::string paraUnit, double value)
        :parameter_name(paraName),parameter_unit(paraUnit),parameter_value(value)
    {

    }
    
};

class ParameterTable
{
public:
    std::vector<para> table;
    int num_of_row;
    std::string table_name;

    ParameterTable(int row=0, std::string name="");
    ParameterTable(std::string filename, std::string start = "", int row = 0);
    ~ParameterTable(){};
    ////////////////////////////////////////////////////////////////////////
    /////////////////////////以下为文件读取操作///////////////////////////////
    //默认表格开始名字为空，即从第一个不为空的行开始读取一直到读到空行
    //默认不知道表格行数，自动判断行数，也可以自己选取读取多少行
    bool readInFile(std::string filename, std::string start = "", int row = 0);


    //////////////////////////以下为输出操作//////////////////////////////
    bool writeToFile(std::string filename, int mode = 0, std::string title = "");
    void show();

    /////////////////////////以下为修改操作//////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    template <class T> bool setParasName(T& str);
    template <class T> bool setParasUnit(T& str);
    template <class T> bool setParasValue(T& arr);
    bool append(std::string paraName, std::string paraUnit, double value);

    /////////////////////////以下为表格查找操作////////////////////////////////
    double searchValue(std::string index);
    std::string searchUnit(std::string index);  

    ///////////////////////// 以下为统计操作////////////////////////////////////
    double sum(int start, int end);
    

private:
    void setUnit2SI();        
};

#endif // !PARAMETERTABLE_H
