#ifndef FUNCTION
#define FUNCTION
#include <iostream>
#include <fstream>
#include <iomanip>
#include <Windows.h>
#include <limits>
#include <stdlib.h>
#include <time.h>
#include <random>
//#include< fstream >

namespace function
{
    //http://www.cplusplus.com/reference/random/
    std::vector<double> randNumber(double start, double end, int number = 100);

    std::string Trim(std::string str);

    /// <summary>
    /// 判断字符串是否为数字
    /// </summary>
    /// <param name="str"></param>
    /// <returns></returns>
    bool isNumber(std::string str);

    /// <summary>
    /// 任意精度打印变量的值
    /// </summary>
    /// <param name="name"></param>
    /// <param name="x"></param>
    /// <param name="precision"></param>
    void print(std::string name, double x, int precision = int(9));

    /// <summary>
    /// 任意数据类型数组的长度
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="array"></param>
    /// <returns></returns>
    template <class T>
    int getArrayLen(T& array);


    std::string getStrUntil(std::ifstream& fin, char deli=',');

    /// <summary>
    /// 字符串转化为double字符
    /// </summary>
    /// <param name="in"></param>
    /// <returns></returns>
    double str2double(std::string in);

    std::string double2string(double x);

    int calStringContainsChar(std::string s, char p);

    std::string to_utf8(const wchar_t* buffer, int len);


    std::string to_utf8(const std::wstring& str);

    std::wstring stringToWstring(const std::string& str);

    //将string转化为utf-8编码以便写入
    std::string to_utf8(const std::string& str);

    //std::string DoubleToString(const double value, unsigned int precisionAfterPoint = 6);
}


#endif // !FUNCTION
