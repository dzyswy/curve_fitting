#include "augmentedmatrix.h"
#include <iostream>

//构造函数1
AugmentedMatrix::AugmentedMatrix()
{
    width = height = 0;
}

//析构函数
AugmentedMatrix::~AugmentedMatrix()
{
    releaseData();
}

//矩阵导入数据
bool AugmentedMatrix::setMatrixData(fraction *valves , int width, int height)
{
    this->width = width;
    this->height = height;

    for(int i=0 ; i<this->height ; i++)
    {
        fraction *x = new fraction[this->width];
        for(int j=0 ; j<this->width ; j++)
            x[j] = valves[i*this->width + j];
        matrix[i] = x;
    }

    return true;
}

// 处理数据
bool AugmentedMatrix::process()
{
    return augmentedMatrixTrans();
}

 // 打印结果
void AugmentedMatrix::print()
{
    std::cout << "\nAugmented Matrix :" << std::endl;
    for( int y = 0 ; y<height ; y++ )
    {
        for( int x = 0 ; x<width ; x++ )
        {
            std::cout << (float)matrix[y][x] << "\t";
        }
        std::cout << (float)augment[y] << std::endl;
    }
}

// 设置增广列数据
bool AugmentedMatrix::setAugmentData( fraction *values , int height)
{
    this->height = height;
    for( int i(0) ; i < this->height ; i++)
        augment[i]  = values[i];
    return true;
}

//增广矩阵消元变换
bool AugmentedMatrix::augmentedMatrixTrans()
{
    if(width != height)
        return false;

    for(int p=0 ; p<width ; p++)
    {
        if(!sortLines(p))// 排序
            return false;
        elimination(p);// 消元
    }

    // 左上角的消元，求根
    for(int x=width-1 ; x>=0 ; x--)
    {
        if(matrix[x][x] == 0.0)
            return false;
        elimination_inv(x);
    }
    identityMatrix();

    return true;
}

// 排序，防止(x,y)坐标处的值为0，p为矩阵的行和列（此时行列相同）
// 返回false，代表方程组无解
bool AugmentedMatrix::sortLines(int p)
{
    for( int y = p ; y<height ; y++ )
    {
        if( matrix[y][p] == 0.0 )
        {
            if( y == height-1 ) //全部为0
                return false;
        }else{
            if( y != p )
                swap(matrix[y] , matrix[p]);
            break;
        }
    }

    return true;
}

// 高斯消元，p为矩阵的行和列（此时行列相同）
void AugmentedMatrix::elimination(int p)
{
    for( int y = p+1 ; y<height ; y++ )
    {
        fraction t = matrix[p][p] / matrix[y][p];
        //matrix[y][p] = 0.0;
        for( int x = p ; x<width ; x++ )
        {
            matrix[y][x] = matrix[y][x]*t - matrix[p][x];
        }
        augment[y] = augment[y]*t - augment[p];
    }
}

// 交换两行
void AugmentedMatrix::swap(fraction *a, fraction *b)
{
    fraction *pF = a;
    a = b;
    a = pF;
}

// 左上角的消元，求根
void AugmentedMatrix::elimination_inv(int x)
{
    for( int y = 0 ; y<x ; y++ )
    {
        if(matrix[y][x] == 0.0)
            continue;

        fraction t = matrix[y][x] / matrix[x][x];
        matrix[y][x] = 0.0;
        augment[y] = augment[y] - augment[x]*t;
    }
}

//化为单位矩阵
void AugmentedMatrix::identityMatrix()
{
    for( int p = 0 ; p<width ; p++ )
    {
        if(matrix[p][p] == 0.0)
            continue;

        augment[p] /= matrix[p][p];
        matrix[p][p] = 1.0;
    }
}

// 释放内存
void AugmentedMatrix::releaseData()
{
     for(int i(0) ; i < height ; i++)
     {
         if(matrix[i] != NULL )
             delete [] matrix[i];
     }

     width = height = 0;
}
