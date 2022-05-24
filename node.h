#pragma once
//对几何数据的存储和预处理
#include <string.h>
#include <mkl.h>
#include <QGraphicsEllipseItem>
#include <set>
#include "global.h"
//节点类
class Node :public QGraphicsEllipseItem {
public:
    
    double X=0;
    double Y=0;
    ULL nodenum = 0;
    double Fx = 0;
    double Fy = 0;
    //父节点序列
    //int* parentElement;
    //ULL parentElement_num = 0;
    std::set<ULL> parentElement;

    Node();
    ~Node() {};
    void set(ULL num, double x, double y);


private:
    
};

