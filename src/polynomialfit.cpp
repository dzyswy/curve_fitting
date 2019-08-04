#include "polynomialfit.h"
#include <math.h>
#include <memory.h>
#include <iostream>

//构造函数1
PolynomialFit::PolynomialFit()
{
    init();
}

// 析构函数
PolynomialFit::~PolynomialFit()
{
    releaseSamples();
    releaseMatrix();
}

/**
 * @brief PolynomialFit::setAttribute 设置多项式的属性
 * @param degree 设置多项式的阶数(degree请不要大于32)
 * @param intercedeEnable 设置是否使能截距
 * @param intercede 多项式曲线的截距
 */
void PolynomialFit::setAttribute(int degree , bool intercedeEnable, double intercede )
{
    this->isIntercedeEnable = intercedeEnable;

    if(this->isIntercedeEnable){
        this->multiterm = degree;
        this->intercede = intercede;
    }else{
        this->multiterm = degree+1;
        this->intercede = 0.0;
    }
    this->degree = degree;

}

//初始化
void PolynomialFit::init()
{
    multiterm = count = degree =0;
    sx =  sy = weight = NULL;

    intercede = 0.0 ;
    isIntercedeEnable = false;

    matrix  = rvector = results = NULL; // 系数矩阵和结果向量
}

/**
 * @brief PolynomialFit::setSample 设置需要拟合曲线的样本，权重等
 * @param x 样本x轴
 * @param y 样本y轴
 * @param count 样本数量
 * @param enableWeight 是否使能权重
 * @param w 样本的权重
 * @return 样本是否可以拟合
 */
bool PolynomialFit::setSample(double *x,
                              double *y,
                              int count ,
                              bool enableWeight ,
                              double *w )
{
    releaseSamples();
    this->count = count;
    this->sx = new fraction[count] , this->sy = new fraction[count] ,\
            this->weight = new fraction[count];

    for(int i(0) ; i < count ; i++)
    {
        this->sx[i] = x[i] ;
        this->sy[i] = y[i] ;
        if(enableWeight)
            this->weight[i] = w[i];
        else
            this->weight[i] = 1.0;
    }

    return true;
}

// 获取计算的结果
float PolynomialFit::getResult(int y)
{
    return (float)results[y];
}

// 数据处理
bool PolynomialFit::process()
{
    // 利用最小二乘法生一次方程组
    if(!leastQuare())
        return false;

    // 通过高斯法解方程组
    AugmentedMatrix am;
    am.setMatrixData( matrix , multiterm , multiterm ) ;
    am.setAugmentData( rvector , multiterm);
    if( !am.process() )
        return false;

    // 解出的结果就是多项式的系数
    for(int i=0;i<multiterm;i++){
        results[i] = am.getResult(i);
    }

    return true;
}

// 打印结果
void PolynomialFit::print()
{
	int y = 0;
    std::cout << "\nPolynomial Fit :" << std::endl;
    for( y = 0 ; y<multiterm ; y++ )
    {
        for( int x = 0 ; x<multiterm ; x++ )
        {
            std::cout << (float)matrix[y*multiterm + x] << "\t";
        }
        std::cout << (float)rvector[y] << std::endl;
    }

    std::cout << "\nresults:" << std::endl;
    for( y = 0 ; y<multiterm ; y++ )
    {
        std::cout << (float)getResult(y) << "\t";
    }
    std::cout << std::endl;
}

// 释放内存
void PolynomialFit::releaseSamples()
{
    if(sx){
        delete [] sx;
        sx = NULL ;
    }
    if(sy){
        delete [] sy;
        sy = NULL ;
    }
    if(weight){
        delete weight;
        weight = NULL;
    }
    count  = 0 ; // 样本数量

}

// 释放矩阵内存
void PolynomialFit::releaseMatrix()
{
    if(matrix ){
        delete [] matrix ;
        matrix = NULL;
    }

    if(rvector){
        delete [] rvector;
        rvector = NULL;
    }

    if(results){
        delete [] results;
        results = NULL;
    }
}

// 利用最小二乘法，将样本转换成求解的方程组
bool PolynomialFit::leastQuare()
{
    int x = 0 , y = 0;
    if(multiterm > count || count ==0 || multiterm == 0)
        return false;

    releaseMatrix();
    matrix = new fraction[multiterm*multiterm];
    rvector = new fraction[multiterm];
    results = new fraction[multiterm];

    //对角矩阵，计算各个行列的值
    for(int i=0 ; i < 2*multiterm-1 ; i++)
    {
        if( i < multiterm){
            x = i;
            y = 0;
        }else{
            x = multiterm - 1;
            y = i - x;
        }

        fraction r = calcMatrixUnit(x, y);
        for(int yi=0;yi<=i;yi++)
        {
            int xi = i-yi;
            if(yi < multiterm && xi < multiterm)
                matrix[yi*multiterm + xi] = r;
        }
    }

    // 计算增广列的值
    for(y=0;y<multiterm;y++)
        fillVectorUnit(y);
    return true;
}

// 计算矩阵某项的累计值
fraction PolynomialFit::calcMatrixUnit(int x, int y)
{
    fraction sum = 0.0;
    double n = 2*degree - x - y  ;
    for(int i=0;i<count;i++)
    {
        sum += pow( (double)sx[i] , n ) * weight[i];
    }

    return sum;
}

// 计算向量某项的累计值
void PolynomialFit::fillVectorUnit(int y)
{
    rvector[y] = 0.0;
    double n = degree - y;
    for(int i=0;i<count;i++)
    {
        rvector[y] += (sy[i]-intercede) * pow( (double)sx[i] , n ) * weight[i];
    }
}

