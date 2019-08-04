#include <stdio.h>
#include "linearfit.h"
#include <iostream>


// 多项式拟合测试
void PolyfitTest()
{
    int degree = 3;
    int count = 5;
    double sx[] = {3.35,100.35,258.23,526.24,935.6};
    double sy[] = {0.5225,10200.8225,67000.1929,277000.0176,697000.56};
    double w[] = {1,0.02,0.01,0.005,0.001};
    Polynomial pf;
    pf.setAttribute(degree , false , 1.0);
    if( pf.setSample(sx,sy,count,false,w) &&
            pf.process() ){
        pf.print();
        std::cout << "\ngoodness: " << pf.getGoodNess() << std::endl;
    }else
        std::cout << "failed" ;
}

// 对数拟合测试
void LogarithmTest()
{
    double sx [] = {1,2,3,4,5,6,7,8,9,10};
    double sy[] = {6,12,24,32,41,50,59,68,77,86};
    Logarithm lf;
    if( lf.setSample(sx,sy,10,false,NULL) &&
            lf.process() ){
        lf.print();
        std::cout << "\ngoodness:" <<  lf.getGoodNess() << std::endl;
    }else
        std::cout << "failed" ;
}

// 指数拟合测试
void ExponentTest()
{
    double sx [] = {1,2,3,4,5};
    double sy[] = {4.190873241,
                   7.636268922,
                   13.91418917,
                   25.35330568,
                   48};
    Exponent ef;
    if( ef.setSample(sx,sy,5,false,NULL) &&
            ef.process() ){
        ef.print();
        std::cout << "\ngoodness:" <<  ef.getGoodNess() << std::endl;
    }else
        std::cout << "failed" ;
}

// 幂函数拟合测试
void PowerTest()
{
    double sx[] = {1,2,3,4,5};
    double sy[] = {1.258,
                   2.615581859,
                   4.013476146,
                   5.438210223,
                   9};
    Power pf;
    if( pf.setSample(sx,sy,5,false,NULL) &&
            pf.process() ){
        pf.print();
        std::cout << "\ngoodness:" << pf.getGoodNess() << std::endl;
    }else
        std::cout << "failed" ;
}

int main( int  , char** )
{
    PolyfitTest();
//    PowerTest();
//    ExponentTest();
//    LogarithmTest();

    return 0;
}

