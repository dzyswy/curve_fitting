//*******************************************************
// filename : augment matrix
// author : lxw
// description : Solving the first order equation group
// date : 20171221
// contact : QQ392134738
//*******************************************************

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <memory.h>
#include <vector>
#include <algorithm>

// 分数的定义--提高计算的精度
/*
class Fraction
{
public:
    double numerator; // 分子
    double denominator; // 分母

    Fraction()
    {
        numerator = 1.0;
        denominator = 1.0;
    }

    Fraction(double numerator , double denominator)
    {
        this->numerator = numerator;
        this->denominator = denominator;
    }

    Fraction(double n)
    {
        numerator = n;
        denominator = 1.0;
    }

    // =
    inline void operator=(double numerator)
    {
        this->numerator = numerator;
        this->denominator = 1.0;
    }

    // =
    inline void operator=(Fraction &a)
    {
        this->numerator = a.numerator;
        this->denominator = a.denominator;
    }


    // -
    friend inline Fraction operator-(Fraction &a , Fraction &b)
    {
        double n = a.numerator*b.denominator - a.denominator*b.numerator;
        double d = a.denominator*b.denominator;
        return Fraction(n , d);
    }

    // +
    friend inline Fraction operator+(Fraction &a , Fraction &b)
    {
        Fraction ret ;
        ret.numerator = a.numerator*b.denominator + a.denominator*b.numerator;
        ret.denominator = a.denominator*b.denominator;
        return ret;
    }

    // -=
    inline void operator-=(Fraction &a)
    {
        numerator = a.denominator*numerator - a.numerator*denominator;
        denominator = a.denominator*denominator;
    }

    // *
    friend inline Fraction operator*(Fraction &a , Fraction &b)
    {
        double n = a.numerator*b.numerator;
        double d = a.denominator*b.denominator;
        return Fraction(n,d);
    }

    // x=
    inline void operator*=(Fraction &a)
    {
        numerator *= a.numerator;
        denominator *= a.denominator;
    }

    // /
    friend inline Fraction operator/(Fraction &a , Fraction &b)
    {
        Fraction ret ;
        ret.numerator = a.numerator*b.denominator;
        ret.denominator = a.denominator*b.numerator;
        return ret;
    }

    // ==
    inline bool operator==(Fraction &a)
    {
        return this->numerator*a.denominator==this->denominator*a.numerator;
    }

    // ==
    inline bool operator==(double a)
    {
        return this->numerator==this->denominator*a;
    }
};
*/

// 使用long double 提高计算的精度
#if !defined (fraction)
typedef long double  fraction;
#endif


// class AugmentedMatrix ： 增广矩阵，利用高斯消元的方法求解
// 一次方程组，计算线性拟合的结果。矩阵的定义-矩阵的坐标如图所示
//-|------------->x
// |
// |
// |
// |
// \/
// y
// 矩阵数据和增广列，如下所示
// |--------matrix-------|--augment--|----result----|
// |--2a + 3b + 11c + 5d =  2 -------|--1 0 0 0|-2--|
// |--1a +  b +  5c + 2d =  1 -------|--0 1 0 0| 0--|
// |--2a +  b +  3c + 2d = -3 -------|--0 0 1 0| 1--|
// |--1a +  b +  3c + 3d = -3 -------|--0 0 0 1|-1--|
class AugmentedMatrix
{
public:
    AugmentedMatrix(); //构造函数1
    ~AugmentedMatrix(); //析构函数

    bool setAugmentData(fraction *values, int height); //矩阵导入数据()
    bool setMatrixData(fraction *valves, int width, int height); // 设置增广列数据()

    inline fraction getResult(int y){return augment[y];}
    bool process(); // 处理数据
    void print(); // 打印结果

private:
    fraction *matrix[32]; // 矩阵数据--定义成分数形式提高精度(最多32组)
    fraction augment[32]; // 增加的一列
    int width , height; // 矩阵的宽高

    void releaseData(); // 释放内存
    bool augmentedMatrixTrans(); // 高斯消元法，求解多元方程组的值
    bool sortLines(int p ); // 排序，
    void elimination( int p ); // 消元
    void swap(fraction *a , fraction *b); // 交换两行
    void elimination_inv( int x ); // 消元
    void identityMatrix(); //化为单位矩阵
};

#endif // _MATRIX_H_
