//******************************************
// filename : polynomial fit
// author : lxw
// description : polynomial fit
// date : 20171221
// contact : QQ392134738
//******************************************


#ifndef _POLYNOMAIAL_FIT_H_
#define _POLYNOMAIAL_FIT_H_

#include "augmentedmatrix.h"

class PolynomialFit
{
public:
    PolynomialFit();
    ~PolynomialFit();

    //设置多项式的属性，如函数的阶数 , 曲线的截距等
    void setAttribute(int degree ,
                      bool intercedeEnable = false,
                      double intercede = 0.0 );

    //设置曲线的样本，权重等
    virtual bool setSample(double *x ,
                           double *y ,
                           int count ,
                           bool enableWeight  ,
                           double *w );

    virtual float getResult(int y); // 获取计算的结果
    virtual bool process();
    void print(); // 打印


protected:
    double intercede; // 截距
    bool isIntercedeEnable; // 截距使能
    int degree; // 阶数
    int multiterm; // 多项式的项数
    fraction *sx , *sy , *weight; // 样本-权重
    int count; // 样本数量

private:
    fraction *matrix , *rvector , *results; // 系数矩阵和结果向量

    bool leastQuare(); // 利用最小二乘法，将样本转换成求解的方程组
    fraction calcMatrixUnit(int x , int y); // 计算矩阵某项的累计值
    void fillVectorUnit(int y); // 计算向量某项的累计值

    void init();
    void releaseSamples(); // 释放内存
    void releaseMatrix();

};

#endif // _POLYNOMAIAL_FIT_H_
