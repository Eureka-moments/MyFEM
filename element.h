#pragma once
//#ifndef ELEMENT_H
//#define ELEMENT_H


//对几何数据的存储和预处理
#include<string.h>
#include<mkl.h>
#include<vector>
#include"node.h"
#include"global.h"
//#include"compute.h"

//extern ULL get_line_index(ULL i, ULL j, ULL n);


class Element {
public:
    //单元编号
    ULL elementnum;
    //从属于单元的节点编号
    ULL node[NODE_NUM_IN_ELEMENT];


    //构造函数
    Element() {
        elementnum = 0;
        memset(node, 0, sizeof(node));
    }
    void set(ULL num, ULL node1, ULL node2, ULL node3, Node* nodeline_input, double* K_input);

    //单元的物理属性
    double E = 0;//弹性模量
    double miu = 0;//泊松比
    double E_scala = 0;
    //double zou=0;//密度


    //单元节点坐标(由位移决定)
    double x[NODE_NUM_IN_ELEMENT * PROBLEM_DIMENSION];//6
    double* y = x + NODE_NUM_IN_ELEMENT;

    //单元节点位移
    double u[NODE_NUM_IN_ELEMENT * PROBLEM_DIMENSION] = { 0 };//6
    double* v = u + NODE_NUM_IN_ELEMENT;

    //单元应变(由位移决定)
    double epsilon[3] = { 0 };//epsilon_x epsilon_y gama_xy

    //单元应力(由位移决定)
    double sigma[3] = { 0 };

    //位移参数(由位移决定)
    double alpha[NODE_NUM_IN_ELEMENT * PROBLEM_DIMENSION] = { 0 };//u=alpha1+alpah2*x+alpha3*y
    double* beta = alpha + NODE_NUM_IN_ELEMENT;//v=beta1+beta2*x+beta3*y_v...

    //单元a、b、c属性(稳定的)
    double a[NODE_NUM_IN_ELEMENT * 3] = { 0 };
    double* b = a + NODE_NUM_IN_ELEMENT;
    double* c = b + NODE_NUM_IN_ELEMENT;

    //单元的面积乘系数(稳定的)
    double S = 0;
    double S_scala = 0;

    //单元K矩阵
    //double K[((PROBLEM_DIMENSION*NODE_NUM_IN_ELEMENT+1)*PROBLEM_DIMENSION*NODE_NUM_IN_ELEMENT)/2]={0};



    //传入参数
    void set_uv(double* Delta);
    void set_E_miu(double E_tmp,double miu_tmp);

    //获取单元内任意一点位移
    double get_u(double x_in, double y_in);
    double get_v(double x_in, double y_in);

private:
    //K矩阵
    double* K=NULL;
	//node_line
    Node* node_line=NULL;
    //xy是否改变
    bool is_xy = false;
    //E、miu是否改变
    bool is_E_miu = false;
    //密度是否改变
    //bool is_zou=false;
    //位移是否改变
    bool is_uv = false;
    //几何参数更新
    void update_abcS();
    //单元应变更新
    void update_epsilon();
    //单元应力更新
    void update_sigma();
    //单元形函数更新
    void update_alpha();
    //E、miu矩阵乘常数
    void update_E_scala();
    //K矩阵更新
    void update_K();

    ULL get_line_index(ULL i, ULL j);
    //检查更新
    bool check_update();

    //形函数计算
    //void compute_shape(double* N, double x, double y);
};



/*#endif*/ // ELEMENT_H
