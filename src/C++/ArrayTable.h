#ifndef ARRAYTABLE_H
#define ARRAYTABLE_H
#include <iostream>
#include <vector>

template <class Type>
class PhysicalVarList
{
public:
    std::string ColName;
    std::string ColUnit;
    std::vector<Type> data;

    PhysicalVarList(int row = 0);

    //运算符<<重载
    friend std::ostream& operator<<(std::ostream& out, PhysicalVarList L);

    bool add(Type num);

    bool Multiply(double num, std::string new_unit = "EmptyUnit");
};

class ArrayTable
{
public:
    unsigned int col, row;
    std::vector<PhysicalVarList<double>> table;
    std::string table_name;

    ArrayTable() {};
    ~ArrayTable() {};
    ArrayTable(int _col, int _row = 0, std::string _tablename = "");
    ArrayTable(std::string _filename, std::string _start = "");

    ////////////////////////输入到表格内/////////////////////////////////
    bool readInFile(std::string _filename, std::string _start="",char delimiter=',');
    //根据文件名和列数读取表格，如果有起始字符则从起始字符Start开始读取，如果没有则从读取到的第一个不为空的行开始读取。
    //应要求自动判断读取的列
    ////////////////////////以上为表格读取数据操作////////////////////////


    ////////////////以下为输出操作//////////////////
    void show();//在命令行打印整个表格
    void show(int column);//展示第column列的内容
    enum class write_to_file_type { append2file, truncate2file };
    bool writeToFile(std::string filename, write_to_file_type mode= write_to_file_type::truncate2file);//输出到文件
    void openWithProgram();
    std::vector<double> selectRowValue(int _row);//读取行row的一条记录
    bool plot(int _coly=1,int _colx=0);
    ///////////////输出和展示操作////////////////////

    ///////////////////////////////////////////////
    ///////////////以下为修改操作////////////////////
    ///////////////////////////////////////////////
    void setUnit2SI();//将单位设置为SI单位
    bool setRowValue(int _row, std::vector<double> _arr);//修改第row列的记录
    bool exchageRowValue(int _row1, int _row2);//交换row1和row2的值
    bool append(std::vector<double> arr);
    bool setParaName(std::vector < std::string> _arr);
    bool setParaUnit(std::vector < std::string> _arr);
    void pop();

    ////////////////清除表格数据操作////////////////////
    void Clear();

    ////////////////////////////////////////////////
    ///////////////以上为修改操作//////////////////////
    ///////////////////////////////////////////////



    ///////////////////查询操作///////////////////////
    ////////////找某一列最大值点对应的行///////////////
    int findMaxDataIndex(int _col);

    int findMinDataIndex(int _col);

    ///////////////找出所有的极大值点对应的行，默认大于前后两个点的值/////////////////
    std::vector<int> findMaximumDataIndex(int _col, int _order = 2);

    int findMaxDataIndex(int _col) const;

    int findMinDataIndex(int _col) const;

    double findMaxData(int _col) const;

    double findMinData(int _col) const;

    ////////////////////////////////////////////////
    ///////////////以下为统计操作/////////////////////
    ////////////////////////////////////////////////
    double sumAll(int column) const;

    double quasierror(int _col, double mean = 0) const;

    void doQuickSort(int _col);//对表格快速排序

    ////////////////////////////////////////////////////////
    ///////////////以下为插值、积分、微分操作////////////////////
    ////////////////////////////////////////////////////////
    double linearInterpolate(double x, int columny = 1,int columnx=0,int type=0);//对第column行线性插值，返回插值结果
                                                                  //type=0,循环插值，type=1不循环,超出范围置为0;
    
    //////////////////////FFT插值函数，具有周期性，但是要求自变量序列是等间隔的，否则不准///////////////
    /////////////////x插值点，x_col差值列x，column y所在列//////////////////////////
    //////////////////////插值效果并不好///////////////////////////////////////////
    double triInterpolate(double _x, int _coly = 1, int _colx = 0);
    double Integrate(int _coly, int _start, int _end,int _colx=0);//column列，对start到end行积分
    void diff(int _coly, int _colx = 0);

private:
    void QuickSort(int _col, ArrayTable point[], int low, int high);
};



#endif // !PARAMETERTABLE_H
