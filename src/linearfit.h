//******************************************
// filename : linear fit
// author : lxw
// description : linear fit
// date : 20171221
// contact : QQ392134738
//******************************************

#ifndef _LINEAR_FIT_H_
#define _LINEAR_FIT_H_


#include "polynomialfit.h"
#include <iostream>
#include <math.h>

/*
 *函数名：GoodNess计算拟合优度R
*/
class GoodNess
{
public:
    GoodNess():goodness(0.0){;}
    ~GoodNess(){;}

    void calcGoodNess(double *x,double *y,int count);
    double getGoodNess();

protected:
    virtual double func(double x) = 0;
    double goodness; // 拟合优度
};


/* Polynomial -- 多项式的拟合y=(Pn*x^n)+(Pn-1*x^n-1)+...+P1*x+P0
 * setAttribute: 设置函数的阶数 , 曲线的截距
 * setSample: 可以设置样本、样本数量、使能权重、每个样本的权重
 * process: 阶数和样本设置后计算系数和拟合优度（必须先设置好多项式的阶数和样本）
*/
class Polynomial : public PolynomialFit , public GoodNess
{
public:
    Polynomial();
    ~Polynomial();

    bool setSample(double *x,
                  double *y,
                  int count ,
                  bool enableWeight ,
                  double *w ); // 设置样本
    bool process(); // 处理，添加拟合优度的计算

protected:
    double *osx , *osy ; // 原始样本内容和数量
    int ocount;
    void copy(double *x, double *y, int count);

    double func(double x); // 方程
};


/* logarithm -- 对数拟合 y = P2*ln(x) + P1
 * setSample: 可以设置样本、样本数量、使能权重、每个样本的权重
 * process: 阶数和样本设置后计算系数和拟合优度（必须先设置好样本）
*/
class Logarithm : public Polynomial
{
public:
    Logarithm();
    bool setSample(double *x,
                  double *y,
                  int count ,
                  bool enableWeight ,
                  double *w ); // 设置样本


    double func(double x);
};


/* Exponent -- 对数拟合 y=P2*exp(P1*x)
 * setSample: 可以设置样本、样本数量、使能权重、每个样本的权重
 * process: 阶数和样本设置后计算系数和拟合优度（必须先设置好样本）
*/
class Exponent : public Polynomial
{
public:
    Exponent();
    bool setSample(double *x,
                  double *y,
                  int count ,
                  bool enableWeight ,
                  double *w ); // 设置样本

    double func(double x);
    float getResult(int y); // 系数从左到右分别对应0、1、2...n
};

/* Power -- 幂函数拟合 y=P2*x^(P1)
 * setSample: 可以设置样本、样本数量、使能权重、每个样本的权重
 * process: 阶数和样本设置后计算系数和拟合优度（必须先设置好样本）
*/
class Power : public Polynomial
{
public:
    Power();
    bool setSample(double *x,
                  double *y,
                  int count ,
                  bool enableWeight ,
                  double *w ); // 设置样本

    double func(double x);
    float getResult(int y); // 系数从左到右分别对应0、1、2...n
};

#endif // _POLYNOMAIAL_FIT_H_
