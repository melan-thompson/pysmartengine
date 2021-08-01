#include<cmath>
#include<iostream>

#ifndef CYLINDERFUNCTIONS
#define CYLINDERFUNCTIONS

/// <summary>
/// 计算平均有效压力
/// </summary>
/// <param name="D"></param>
/// <param name="cm"></param>
/// <param name="pme"></param>
/// <param name="stroke_num"></param>
/// <returns></returns>
double FMEP(double D,double cm,double pme,int stroke_num=4)
{
    if(stroke_num==4)
    {
        return pow(D, -0.1778) * (0.0855 * cm + 0.0789 * pme / 98066.5 - 0.214) * 98066.5;
    }
    else if(stroke_num==2)
    {
        return pow(D, -0.1988) * (0.0829 * cm + 0.08295 * (pme / 98066.5) - 0.3) * 98066.5;
    }
    else
    {
        using namespace std;
        std::cerr << "number of stroke error! only 2 or 4 is allowed!\n";
    }

}

#endif // !CYLINDERFUNCTIONS