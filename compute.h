#pragma once
//#ifndef COMPUTE_H
//#define COMPUTE_H

//计算函数头文件——所有矩阵变化和解方程
#include "mkl.h"
#include "element.h"
#include "node.h"
#include "global.h"
#include <string.h>
#include <iostream>


//根据方程解更新原向量
bool update_original_delta(ULL element_num, ULL constraint_num, ULL* constraint_level, double* F_resize, double* Delta, Element* element_line);
bool update_resize_F(ULL constraint_num, ULL* constraint_level, double* F, double* F_resize);
//根据约束创建全变量K矩阵
bool update_resize_K(ULL constraint_num, ULL* constraint_line, ULL* constraint_level, double* K, double* K_resize);

bool compute_K_Delta_F(ULL element_num, ULL node_num, ULL constraint_num, ULL* constraint_level, double* K_resize, double* F_resize, double* Delta, Element* element_line);

//简单的向量算法
double vector_dot_product(double *a, double *b, int n);
void add_vector(double *a, double *b, double *c, int n);
void subtract_vector(double *a, double *b, double *c, int n);
void multiply_scalar(double s, double* a, int n);
double vecor_length(double *a,double *b,int n);

//任务目标实现
//1.计算外力分布向量
// （坐标与节点/单元的对应问题。）
//2.计算单元刚度矩阵
//3.合并计算总刚矩阵
//4.求解总刚的位移-外力大矩阵线性方程
//5.求得位移后回代得到其他一些力学分析信息

/*#endif */// COMPUTE_H
