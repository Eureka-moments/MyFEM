#pragma once
#ifndef TEST_H
#define TEST_H

#include <QtWidgets/QMainWindow>
#include "ui_Test.h"
#include <QFileDialog>
#include "compute.h"
#include "global.h"
#include "GraphicsTriItem.h"
#include "GraphicsCirItem.h"
#include <QGraphicsPolygonItem>
#include <QMessageBox>
#include <QButtonGroup>
#include<QMouseEvent>
#include <mkl.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <set>
#include <iostream>
#include <fstream>



class Test : public QMainWindow
{
    Q_OBJECT

public:
    Test(QWidget *parent = Q_NULLPTR);
	~Test();
	//void find_element(QGraphicsItem* item);
	bool set_zou(double zou_input) {
		zou= zou_input;
		return true;
	};

	QGraphicsScene* scene;
	QGraphicsScene* scene_behind;
	QGraphicsScene* scene_new;
	
	GraphicsTriItem* tri_item;
	GraphicsTriItem* tri_item_new;
	
	//QGraphicsEllipseItem* cons_item;
	//QGraphicsLineItem* force_item;

	ULL element_num = 1;//单元总数(考虑写成全局)
	ULL node_num = 1;//节点总数(考虑写成全局)
	ULL constraint_num = 1;//约束总数(考虑写成全局)
	//输入节点矩阵（节点编号，x，y坐标）?不一定需要这样申请成线性内存
	Node* node_line;
	Element* element_line;
	//约束序号数组(从0开始)
	ULL* constraint_line;
	//由约束得出的层数分割
	ULL* constraint_level;
	//位移向量(uv) 一定需要这样申请成线性内存（大矩阵乘法）
	double* Delta;
	//resize后的位移向量
	//double* Delta_resize;
	//节点作用力?不一定需要这样申请成线性内存
	double* F;
	//resize后的F向量
	double* F_resize;
	//总刚矩阵存储 一定需要这样申请成线性内存（大矩阵乘法）
	double* K;
	//修正后的刚度矩阵
	double* K_resize;
	
	
public slots:
	void get_QPointF(QPointF point);
	void close_click();
	
private slots:

	bool choose_model_file();//节点、单元、材料
	bool choose_boundary_file();//约束、外力
	
	void confirm_click();
	void clean_click();
	void color_change();
	void get_two_number();
	void check_QPointF_information(double x,double y,QList<QGraphicsItem*> element_items);
	void check_point_information(double x, double y, QList<QGraphicsItem*> element_items, QList<QGraphicsItem*> node_items);
	
private:
	double zou = 0.007;//密度
	bool if_input_model = false;
	bool if_input_boundary = false;
	QButtonGroup* block1;
	
	bool Input_Model(char* model_location);
	bool Input_Boundary_file(char* boundary_file);
	bool draw_element();
	bool draw_node();
	bool redraw();
	bool select_redraw();
	
	////对所有单元的应力应变信息进行统计，在正负域分别取到最大最小和其范围
	/*double max_poistive_epsilon[12] = { 0 };
	double measure_len[12] = { 0 };
	double* max_negtive_epsilon = max_poistive_epsilon + 3;
	double* max_poistive_sigma = max_negtive_epsilon + 3;
	double* max_negtive_sigma = max_poistive_sigma + 3;*/
	double positive_mean[6] = { 0 };
	double negtive_mean[6] = { 0 };
	double positive_std[6] = { 0 };
	double negtive_std[6] = { 0 };

	//人为改变指标
	bool is_e_s = 0;//选择指标为应力还是应变True为应变
	int select_kind = 0;//指标对应序号
	unsigned int constrast = 1;//颜色强调强度
	bool boundary_color = false;//边界是否也着色
	int toumingdu = 255;//透明度
	double concentrate_rate = 0.0001;
	
private:
    Ui::TestClass ui;
};


#endif