#include "Functions.h"
#include <string>
#include <sstream>

std::vector<double> function::randNumber(double start, double end, int number)
{
    std::random_device generator;
    std::uniform_real_distribution<double> distribution(start, end);

    std::vector<double> result;
    for (int i = 0; i < number; i++)
    {
        double dice_roll = distribution(generator);
        result.push_back(dice_roll);
    }

    return result;
}


/// <summary>
/// 删除首末的特定字符
/// </summary>
/// <param name="str"></param>
/// <returns></returns>
std::string function::Trim(std::string str)
{
    std::string blanks(" \n,");
    str.erase(0, str.find_first_not_of(blanks));
    str.erase(str.find_last_not_of(blanks) + 1);
    return str;
}


bool function::isNumber(std::string str)
{
    std::stringstream temp(str);
    double t;
    char p;
    if (!(temp>>t))
    {
        return false;
    }
    if (temp>>p)
    {
        return false;
    }
    else
    {
        return true;
    }

}


void function::print(std::string name, double x, int precision)
{
    using namespace std;
    cout << setw(4) << name << " = " << setw(int(precision + int(1))) << setiosflags(ios::fixed) << setprecision(precision) << x << ";" << endl;
    cout << endl;
}


template <class T>
int function::getArrayLen(T& array)
{

    return (sizeof(array) / sizeof(array[0]));

};

std::string function::getStrUntil(std::ifstream& fin, char deli)
{
    long pos = fin.tellg();
    std::string result;
    std::getline(fin, result, deli);
    //if (result.find("\n") != std::string::npos)//含有换行符，重新读取
    //{
    //    fin.seekg(pos, std::ios::beg);
    //    std::getline(fin, result,'\n');
    //}
    return Trim(result);

    
};

double function::str2double(std::string in)
{
    if (!isNumber(in))
    {
        std::cerr << "The string:" << in << " is not a number!!!";
        std::cerr << "Will return -inf.\n";
        system("pause");
        return -INFINITY;
    }
    double result;
    std::stringstream ss;
    ss << Trim(in);
    ss >> result;
    return result;
}

std::string function::double2string(double x)
{
    std::stringstream ss;
    ss << x;
    std::string result;
    ss >> result;
    return result;
}



std::string function::to_utf8(const wchar_t* buffer, int len)
{
    using namespace std;
    int nChars = ::WideCharToMultiByte(
        CP_UTF8,
        0,
        buffer,
        len,
        NULL,
        0,
        NULL,
        NULL);
    if (nChars == 0)return"";

    string newbuffer;
    newbuffer.resize(nChars);
    ::WideCharToMultiByte(
        CP_UTF8,
        0,
        buffer,
        len,
        const_cast<char*>(newbuffer.c_str()),
        nChars,
        NULL,
        NULL);

    return newbuffer;
}


std::string function::to_utf8(const std::wstring& str)
{
    return to_utf8(str.c_str(), (int)str.size());
}

std::wstring function::stringToWstring(const std::string& str)
{
    LPCSTR pszSrc = str.c_str();
    int nLen = MultiByteToWideChar(CP_ACP, 0, pszSrc, -1, NULL, 0);
    if (nLen == 0)
        return std::wstring(L"");

    wchar_t* pwszDst = new wchar_t[nLen];
    if (!pwszDst)
        return std::wstring(L"");

    MultiByteToWideChar(CP_ACP, 0, pszSrc, -1, pwszDst, nLen);
    std::wstring wstr(pwszDst);
    delete[] pwszDst;
    pwszDst = NULL;

    return wstr;
}

//将string转化为utf-8编码以便写入
std::string function::to_utf8(const std::string& str)
{
    std::wstring ws = stringToWstring(str);
    return to_utf8(ws);
}

int function::calStringContainsChar(std::string s, char p)
{
    int result = 0;
    s = Trim(s);
    for (int i = 0; i < s.size(); i++)
    {
        if (s[i] == p)result++;
    }
    return result;
}

//std::string function::DoubleToString(const double value, unsigned int precisionAfterPoint = 6)
//{
//    std::ostringstream out;
//    // 清除默认精度
//    out.precision(std::numeric_limits<double>::digits10);
//    out << value;
//
//    std::string res = std::move(out.str());
//    auto pos = res.find('.');
//    if (pos == std::string::npos)
//        return res;
//
//    auto splitLen = pos + 1 + precisionAfterPoint;
//    if (res.size() <= splitLen)
//        return res;
//
//    return res.substr(0, splitLen);
//}