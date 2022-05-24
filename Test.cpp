#include "Test.h"


Test::Test(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
	this->setWindowTitle(tr("三节点有限单元法平面求解软件v1.0"));

	scene = new QGraphicsScene(this);
	scene_behind = new QGraphicsScene(this);
	//ui.graphicsView->setScene(scene_behind);
	ui.graphicsView->setScene(scene);
	
	block1 = new QButtonGroup(this);
	block1->addButton(ui.radioButton, 0); 
	block1->addButton(ui.radioButton_2, 1);
	block1->addButton(ui.radioButton_3, 2);
	block1->addButton(ui.radioButton_4, 3);
	block1->addButton(ui.radioButton_5, 4);
	block1->addButton(ui.radioButton_6, 5);
	block1->setExclusive(true);
	
	ui.groupBox->hide();
	
	connect(ui.pushButton_3, &QPushButton::clicked, this, &Test::get_two_number);
	connect(block1, &QButtonGroup::buttonClicked, this, &Test::color_change);
	connect(ui.toolButton_5, &QToolButton::clicked, this, &Test::choose_model_file);
	connect(ui.horizontalSlider, &QSlider::valueChanged, this, &Test::color_change);
	connect(ui.horizontalSlider_2, &QSlider::valueChanged, this, &Test::color_change);
	connect(ui.checkBox, &QCheckBox::stateChanged, this, &Test::color_change);
	connect(ui.toolButton_9, &QToolButton::clicked, this, &Test::choose_boundary_file);
	connect(ui.pushButton, &QPushButton::clicked, this, &Test::confirm_click);
	connect(ui.pushButton_2, &QPushButton::clicked, this, &Test::clean_click);
}



bool Test::choose_model_file() {
	if (if_input_model) {
		QMessageBox::warning(this, tr("警告"), tr("已经输入"));
		return false;
	}
	QString file_name = QFileDialog::getOpenFileName(this, tr("选择文件"), ".", tr("文本文件(*.txt)"));
	ui.lineEdit_5->setText(file_name);//修改对应交互输入槽

	if (file_name.isEmpty()) {
		QMessageBox::warning(this, tr("警告"), tr("请输入文件名"));
		if_input_model = false;
		return false;
	}
	QFile file(file_name);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		QMessageBox::warning(this, tr("读取文件"), tr("无法读取文件 %1:\n%2.").arg(file_name).arg(file.errorString()));
		if_input_model = false;
		return false;
	}
	file.close();
	char* model_location;
	QByteArray ba = file_name.toLatin1();
	model_location = ba.data();
	clock_t start, finish;
	double total;
	start = clock();
	if (Input_Model(model_location)) {
		tri_item = new GraphicsTriItem[element_num];
		ui.textBrowser->insertPlainText("已经输入节点、单元、材料文件！\n");
		if (draw_node()&&draw_element()) {
			finish = clock();
			total = (double)(finish - start) / CLOCKS_PER_SEC;
			QString tmp;
			ui.textBrowser->insertPlainText("单元、节点绘制成功！\n");
			ui.textBrowser->insertPlainText(tmp.asprintf("CPU占用的总时间：%fs\n", total));
			
		}
		else {
			ui.textBrowser->insertPlainText("单元、节点绘制失败！\n");
			if_input_model = false;
			return false;
		}
	}
	else {
		ui.textBrowser->insertPlainText("输入模型文件失败！\n");
		if_input_model = false;
		return false;
	}
	if_input_model = true;
	return true;
}
bool Test::choose_boundary_file() {
	if (if_input_boundary) {
		QMessageBox::warning(this, tr("警告"), tr("已经输入"));
		return false;
	}
	QString file_name = QFileDialog::getOpenFileName(this, tr("选择文件"), ".", tr("文本文件(*.txt)"));
	ui.lineEdit_9->setText(file_name);//修改对应交互输入槽

	if (file_name.isEmpty()) {
		QMessageBox::warning(this, tr("警告"), tr("请输入文件名"));
		if_input_boundary = false;
		return false;
	}
	QFile file(file_name);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		QMessageBox::warning(this, tr("读取文件"), tr("无法读取文件 %1:\n%2.").arg(file_name).arg(file.errorString()));
		if_input_boundary = false;
		return false;
	}
	file.close();
	char* bound_location;
	QByteArray ba = file_name.toLatin1();
	bound_location = ba.data();
	clock_t start, finish;
	double total;
	start = clock();
	if (Input_Boundary_file(bound_location)) {
		finish = clock();
		total = (double)(finish - start) / CLOCKS_PER_SEC;
		QString tmp;
		ui.textBrowser->insertPlainText("输入约束、外力作用成功！\n");
		ui.textBrowser->insertPlainText(tmp.asprintf("CPU占用的总时间：%fs\n", total));
	}
	else {
		ui.textBrowser->insertPlainText("输入约束、外力作用失败！\n");
		if_input_boundary = false;
		return false;
	}
	if_input_boundary = true;
	return true;
}
bool Test::Input_Boundary_file(char* boundary_file) {
	
	//input_constraint_file
	
	std::ifstream inputfile;
	inputfile.open(boundary_file);
	ULL row_num = 0;
	inputfile >> constraint_num >> row_num;
	constraint_line = new ULL[constraint_num + 2];
	constraint_level = new ULL[constraint_num + 2]; //层内行数
	F_resize = (double*)calloc(2 * node_num - constraint_num, SIZE_DOUBLE);
	K_resize = (double*)calloc(((2 * node_num - constraint_num + 1) * (2 * node_num - constraint_num)) / 2, SIZE_DOUBLE);
	//cons_item = new QGraphicsEllipseItem[row_num];

	ULL tmp_constraint_num = 1;
	//初始化约束
	//约束数据规定：x坐标、y坐标、x,y方向约束（0/1））
	double x_tmp, y_tmp;
	int cons_dir;
	//QTransform transform;
	QList<QGraphicsItem*> items;
	QPen pen = QPen(Qt::red);
	pen.setColor(Qt::red);
	for (ULL i = 0; i < row_num; i++) {
		inputfile >> x_tmp >> y_tmp >> cons_dir;
		//点触法
		items = scene_behind->items(QPointF(x_tmp * EXPEND_COEF, -y_tmp * EXPEND_COEF));
		if (items.isEmpty()) {
			QString tmp;
			ui.textBrowser->insertPlainText(tmp.asprintf("约束文件中的点(%0.2f,%0.2f)不在图形中！\n", x_tmp, y_tmp));
			return false;
		}
		ULL no_num = 0;
		double distance = 1000;
		for (int s = 0; s < items.size(); s++) {
			Node* select_item = (Node*)items[s];
			select_item->setBrush(Qt::red);
			double local_distance = fabs(select_item->X - x_tmp) + fabs(select_item->Y - y_tmp);
			if (local_distance >= distance)
				continue;
			else {
				distance = local_distance;
				no_num = select_item->nodenum;
			}
		}


		if (no_num == 0) {
			QString tmp;
			ui.textBrowser->insertPlainText(tmp.asprintf("约束点(%0.2f,%0.2f)未作用，出错！\n", x_tmp, y_tmp));
			return false;
		}
		else {
			if (cons_dir == 0) {
				//x方向约束
				constraint_line[tmp_constraint_num] = 2 * (no_num - 1);
				tmp_constraint_num++;
				QString tmp;
				ui.textBrowser->insertPlainText(tmp.asprintf("约束点(%0.2f,%0.2f)x方向约束已作用\n", x_tmp, y_tmp));
			}
			else if (cons_dir == 1) {
				//y方向约束
				constraint_line[tmp_constraint_num] = 2 * (no_num - 1) + 1;
				tmp_constraint_num++;
				QString tmp;
				ui.textBrowser->insertPlainText(tmp.asprintf("约束点(%0.2f,%0.2f)y方向约束已作用\n", x_tmp, y_tmp));
			}
			else if (cons_dir == 2) {
				//xy方向约束
				constraint_line[tmp_constraint_num] = 2 * (no_num - 1);
				tmp_constraint_num++;
				constraint_line[tmp_constraint_num] = 2 * (no_num - 1) + 1;
				tmp_constraint_num++;
				QString tmp;
				ui.textBrowser->insertPlainText(tmp.asprintf("约束点(%0.2f,%0.2f)xy方向约束已作用\n", x_tmp, y_tmp));
			}
			else {
				QString tmp;
				ui.textBrowser->insertPlainText(tmp.asprintf("约束点(%0.2f,%0.2f)未作用，出错！\n", x_tmp, y_tmp));
				return false;
			}
		}
	}
	//约束从小到大排序
	ULL* constraint_line_tmp = new ULL[constraint_num];
	memcpy(constraint_line_tmp, constraint_line + 1, constraint_num * sizeof(ULL));
	std::sort(constraint_line_tmp, constraint_line_tmp + constraint_num);
	memcpy(constraint_line + 1, constraint_line_tmp, constraint_num * sizeof(ULL));
	delete[] constraint_line_tmp;

	constraint_line[0] = -1;
	constraint_line[constraint_num + 1] = 2 * node_num;

	for (ULL i = 1; i <= constraint_num + 1; i++)
		constraint_level[i] = constraint_line[i] - constraint_line[i - 1] - 1;

	
	//更新K矩阵
	if (update_resize_K(constraint_num, constraint_line, constraint_level, K, K_resize)) {
		ui.textBrowser->insertPlainText("K矩阵已重组\n");
	}
	else {
		ui.textBrowser->insertPlainText("K矩阵未重组\n");
		return false;
	}

	//input_force_file
	ULL F_num = 0;//作用数量
	double f_vec[2] = { 0 };//作用大小向量
	double start_point[2] = { 0 };//作用起始点
	double end_point[2] = { 0 };//作用中止点
	QPointF start_QPointF, end_QPointF;
	QPainterPath path;//作用路径
	int f_kind;//作用类型


	//外力数据规定：xy方向f大小、作用起点xy坐标、作用终点xy坐标、场力1或表面力0
	QString tmp;
	inputfile >> F_num;
	for (ULL i = 0; i < F_num; i++) {
		//输入作用信息
		inputfile >> f_vec[0] >> f_vec[1] >> start_point[0] >> start_point[1] >> end_point[0] >> end_point[1] >> f_kind;
		//QPointF作用起始点，作用中止点
		start_QPointF = QPointF(EXPEND_COEF * start_point[0], -EXPEND_COEF * start_point[1]);
		end_QPointF = QPointF(EXPEND_COEF * end_point[0], -EXPEND_COEF * end_point[1]);

		double path_vec[2] = { end_point[1] - start_point[1],end_point[0] - start_point[0] };//路径向量
		double path_offset = vector_dot_product(start_point, path_vec, 2);//计算直线交点的常数项
		double normal_vec[2] = { path_vec[0],-path_vec[1] };//法向量
		path.moveTo(start_QPointF);
		path.lineTo(end_QPointF);//规划处路径

		//根据作用类型修改F向量
		if (f_kind == 1) {
			//场力
			double scala_tmp;
			for (ULL j = 1; j <= element_num; j++) {
				scala_tmp = zou * t * element_line[j].S / 3;
				for (int k = 0; k < 3; k++)
					cblas_daxpy(2, scala_tmp, f_vec, 1, F + (element_line[j].node[k] - 1) * 2, 1);
				/*F[(element_line[j].node[k]-1) * 2] += (1 / 3) * zou * f_x * t * element_line[j].S;
				F[1 + (element_line[j].node[k]-1) * 2] += (1 / 3) * zou * f_y * t * element_line[j].S;*/
			}
			ui.textBrowser->insertPlainText(tmp.asprintf("场力（%0.2f,%0.2f）已经作用！\n", f_vec[0], f_vec[1]));
		}
		else {//表面力
			QList<QGraphicsItem*> element_items;//碰撞到的单元
			QList<QGraphicsItem*> node_items;//碰撞到的节点
			if (fabs(start_point[0] - end_point[0]) < error && fabs(start_point[1] - end_point[1]) < error) {//集中力
				ui.textBrowser->insertPlainText(tmp.asprintf("第%d作用力为点集中力，作用点(%.2f,%.2f)\n", i, start_point[0], start_point[1]));
				element_items = scene->items(start_QPointF);
				node_items = scene_behind->items(start_QPointF);
				//寻找最近节点距离
				double min_distance = 1000000;
				int min_node_index = 0;
				double tmp_distance = 0;
				if (!node_items.isEmpty()) {
					for (int j = 0; j < node_items.size(); j++) {
						Node* select_node_item = (Node*)node_items[j];
						tmp_distance = pow(select_node_item->X - start_point[0], 2) + pow(select_node_item->Y - start_point[1], 2);
						if (min_distance > tmp_distance) {
							min_distance = tmp_distance;
							min_node_index = j;
						}
					}
				}
				//判断作用位置
				if (element_items.isEmpty() && node_items.isEmpty()) {//如果作用点不在单元内，作用点也不在节点上
					ui.textBrowser->insertPlainText(tmp.asprintf("作用点(%.2f,%.2f)未在模型中，作用力未作用！\n", start_point[0], start_point[1]));
					continue;
				}
				else if (node_items.isEmpty() || min_distance > 1) {//如果作用点在单元内,作用点不在节点上或者离节点很远
					GraphicsTriItem* select_element_item = (GraphicsTriItem*)element_items[0];
					ui.textBrowser->insertPlainText(tmp.asprintf("作用力点(%.2f,%.2f)大小(%.2f,%.2f)已经作用在单元上\n", start_point[0], start_point[1], f_vec[0], f_vec[1]));
					select_element_item->setBrush(Qt::black);
					ULL ele_num = select_element_item->element_id;
					double shape_vec[3] = { 1,start_point[0],start_point[1] };
					double N[3] = { 0 };
					cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, element_line[ele_num].S_scala, element_line[ele_num].a, 3, shape_vec, 1, 0, N, 1);
					for (int j = 0; j < 3; j++)
						cblas_daxpy(2, N[j], f_vec, 1, F + (element_line[ele_num].node[j] - 1) * 2, 1);
					/*F[(element_line[ele_num].node[j] - 1) * 2] += f_x * N[j];
					F[1 + (element_line[ele_num].node[j] - 1) * 2] += f_y * N[j];*/
				}
				else {//作用点在节点上
					Node* select_node_item = (Node*)node_items[min_node_index];
					ui.textBrowser->insertPlainText(tmp.asprintf("作用力点(%.2f,%.2f)大小(%.2f,%.2f)已经作用在节点(%.2f,%.2f)上\n", start_point[0], start_point[1], f_vec[0], f_vec[1], select_node_item->X, select_node_item->Y));
					select_node_item->setBrush(Qt::black);
					cblas_daxpy(2, 1, f_vec, 1, F + (select_node_item->nodenum - 1) * 2, 1);
				}

			}
			else {//分布力
				double shape_point_vec[3] = { 1,0,0 };//由形函数得到的点的向量[1,x,y]
				double* target_point = shape_point_vec + 1;//目标简化的集中力作用点
				double target_f[2] = { 0 };
				double start_tmp[4] = { 0 };
				double* end_tmp = start_tmp + 2;//选中单元内作用路径起始点
				element_items = scene->items(path);
				node_items = scene_behind->items(path);
				if (element_items.isEmpty() && node_items.isEmpty()) {//如果作用点不在单元内，作用点也不在节点上
					ui.textBrowser->insertPlainText(tmp.asprintf("作用%d未在模型中，作用力未作用！\n", i));
					continue;
				}
				else if (element_items.isEmpty() && node_items.size() >= 2) {
					/*只碰撞到了点：
					1.整个路径在模型外，只是擦过一个点
					2.在模型边界作用，但路径短，只碰到一个点
					3.在模型边界作用，路径长，碰到很多点
					我们假设单元足够密集，作用路径足够长在边界至少也会经过两个节点,即其他情况我们认为没有作用且经过节点之外两端不考虑
					任务：寻找两相邻碰撞节点所在的公共单元*/
					Node* pre_select_node_item, * next_select_node_item;
					for (ULL j = 1; j < node_items.size(); j++) {
						pre_select_node_item = (Node*)node_items[j - 1];
						next_select_node_item = (Node*)node_items[j];
						start_tmp[0] = pre_select_node_item->X;
						start_tmp[1] = pre_select_node_item->Y;
						end_tmp[0] = next_select_node_item->X;
						end_tmp[1] = next_select_node_item->Y;
						std::set<ULL> out;
						std::set_intersection(pre_select_node_item->parentElement.begin(), pre_select_node_item->parentElement.end(), next_select_node_item->parentElement.begin(), next_select_node_item->parentElement.end(), inserter(out, out.begin()));
						if (out.size() == 0) {//没有交集单元
							ui.textBrowser->insertPlainText(tmp.asprintf("作用%d作用第%d阶段未在模型单元上，作用力未作用！\n", i, j));
							continue;
						}
						//至少一个单元
						double scala_tmp = vecor_length(start_tmp, end_tmp, 2);
						add_vector(start_tmp, end_tmp, target_point, 2);
						multiply_scalar(0.5, target_point, 2);
						cblas_daxpy(2, scala_tmp, f_vec, 1, target_f, 1);
						/*for(int k=0;k<2;k++)
							target_f[k] = f_vec[k] * scala_tmp;*/
						cblas_daxpy(2, 1, target_f, 1, F + (pre_select_node_item->nodenum - 1) * 2, 1);
						cblas_daxpy(2, 1, target_f, 1, F + (next_select_node_item->nodenum - 1) * 2, 1);
					}

				}
				else if (!element_items.isEmpty()) {//碰到了单元,在模型内

					GraphicsTriItem* select_element_item;//当前选中的单元
					Element* select_element;
					//ULL ele_num;//当前选中的单元号

					for (ULL j = 0; j < element_items.size(); j++) {
						int on_line[9] = { 0 };
						int* above_line = on_line + 3;
						int* under_line = above_line + 3;
						int on_line_num = 0;
						int above_line_num = 0;
						int under_line_num = 0;
						memset(target_point, 0, 2 * sizeof(double));//目标简化的集中力作用点
						memset(target_f, 0, 2 * sizeof(double));
						//memset(start_tmp, 0, 4 * sizeof(double));

						select_element_item = (GraphicsTriItem*)element_items[j];
						select_element_item->setBrush(Qt::black);
						//ele_num = select_element_item->element_id;
						select_element = element_line + select_element_item->element_id;
						//得到每个经过单元中的每一段路径的起始点
						for (int k = 0; k < 3; k++) {
							double node_vec[2] = { select_element->x[k] - start_point[0],select_element->y[k] - start_point[1] };//点
							double dot_product = vector_dot_product(normal_vec, node_vec, 2);
							if (fabs(dot_product) < error) {//近似认为穿过该点
								on_line[on_line_num] = k;
								start_tmp[2 * on_line_num] = select_element->x[k];
								start_tmp[2 * on_line_num + 1] = select_element->y[k];
								on_line_num++;
							}
							else if (dot_product > 0) {//点在直线上
								above_line[above_line_num] = k;
								above_line_num++;
							}
							else {
								under_line[under_line_num] = k;
								under_line_num++;
							}
						}

						//计算不同相交情况下的直线向量与求解方程右常数
						double line_vec[4] = { 0 };
						double line_offset[2] = { 0 };
						int cross_point_num = 0;
						
						if (on_line_num == 2) {}//两点在线上
						else {
							if (on_line_num == 1) {//一点在直线上
								if (above_line_num == 1)//确保路径不是边缘擦过，而是过第三边
									cross_point_num = 1;//计算两线交点作为end_tmp
								else//边缘擦过
									continue;
							}
							else if (on_line_num == 0) {//穿过两边
								cross_point_num = 2;
								//计算三线两交点作为start_tmp/end_tmp
							}
							double line_start_point[2] = { select_element->x[above_line[3 * (above_line_num - 1)]] ,select_element->y[above_line[3 * (above_line_num - 1)]] };

							for (int k = 0; k < cross_point_num; k++) {
								double line_tmp_vec[2] = { line_start_point[1] - select_element->y[above_line[3 * (2 - above_line_num) + k]],line_start_point[0] - select_element->x[above_line[3 * (2 - above_line_num) + k]] };
								memcpy(line_vec + 2 * k, line_tmp_vec, 2 * sizeof(double));
								line_offset[k] = vector_dot_product(line_start_point, line_tmp_vec, 2);
							}

							//计算起始点
							double A[4] = { 0 };
							double b[2] = { 0 };
							for (int k = 0; k < cross_point_num; k++) {
								memcpy(A, path_vec, 2 * sizeof(double));
								memcpy(A + 2, line_vec + 2 * k, 2 * sizeof(double));
								b[0] = path_offset;
								b[1] = line_offset[k];
								int ipiv[2] = { 0 };
								LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 2, 2, A, 2, ipiv);
								LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', 2, 1, A, 2, ipiv, b, 1);
								memcpy(start_tmp + 2 * (on_line_num + k), b, 2 * sizeof(double));
							}

							if (j == 0 || j == element_items.size() - 1) {
								double* stay_point, *another_point;
								
								if (j == 0) {
									stay_point = start_point;
									another_point = end_point;
								}else {
									stay_point = end_point;
									another_point = start_point;
								}
								
								double* move_point;
								
								double start_tmp_len = vecor_length(another_point, start_tmp, 2);
								double end_tmp_len = vecor_length(another_point, end_tmp, 2);
								
								if (start_tmp_len>end_tmp_len)
									move_point = start_tmp;
								else
									move_point = end_tmp;
								
								double far_tmp[2] = { 0 };
								double line_tmp[2] = { 0 };

								subtract_vector(stay_point, another_point, line_tmp, 2);
								
								subtract_vector(move_point, stay_point, far_tmp, 2);
								double eigen_value = vector_dot_product(far_tmp, line_tmp, 2);
								
								if (eigen_value > 0)
									memcpy(move_point, stay_point, 2 * sizeof(double));
								
							}
						}
						double scala_tmp = vecor_length(start_tmp, end_tmp, 2);
						add_vector(start_tmp, end_tmp, target_point, 2);
						multiply_scalar(0.5, target_point, 2);
						cblas_daxpy(2, scala_tmp, f_vec, 1, target_f, 1);
							
						//加载到节点上
						double N[3] = { 0 };
						cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, select_element->S_scala, select_element->a, 3, shape_point_vec, 1, 0, N, 1);
						for (int k = 0; k < 3; k++) {
							cblas_daxpy(2, N[k], target_f, 1, F + (select_element->node[k] - 1) * 2, 1);
						}
					}
				}
			}
		}
	}
	if (update_resize_F(constraint_num, constraint_level, F, F_resize)) {
		ui.textBrowser->insertPlainText("F向量已裁剪\n");
	}
	else {
		ui.textBrowser->insertPlainText("F向量未裁剪\n");
		return false;
	}
	
	inputfile.close();
	return true;
}
bool Test::Input_Model(char* model_location) {
	//input_node
	ULL temp_num;
	double temp_x, temp_y;
	std::ifstream inputfile;
	inputfile.open(model_location);
	inputfile >> node_num;

	node_line = new Node[node_num + 1];
	Delta = (double*)calloc(2 * node_num, SIZE_DOUBLE);
	F = (double*)calloc(2 * node_num, SIZE_DOUBLE);
	K = (double*)calloc(((2 * node_num + 1) * node_num), SIZE_DOUBLE);

	//初始化点
	for (ULL i = 1; i <= node_num; i++) {
		inputfile >> temp_num >> temp_x >> temp_y;
		node_line[i].set(temp_num, temp_x, temp_y);
	}



	//input_element
	ULL temp_node[3];
	/*ULL temp_num;*/

	inputfile >> element_num;
	element_line = new Element[element_num + 1];

	//初始化单元
	for (ULL i = 1; i <= element_num; i++) {
		inputfile >> temp_num >> temp_node[0] >> temp_node[1] >> temp_node[2];
		double x_1 = node_line[temp_node[1]].X - node_line[temp_node[0]].X;
		double x_2 = node_line[temp_node[2]].X - node_line[temp_node[0]].X;
		double y_1 = node_line[temp_node[1]].Y - node_line[temp_node[0]].Y;
		double y_2 = node_line[temp_node[2]].Y - node_line[temp_node[0]].Y;
		if (x_1 * y_2 - x_2 * y_1 < 0) {
			element_line[i].set(temp_num, temp_node[0], temp_node[2], temp_node[1],node_line,K);
		}
		else {
			element_line[i].set(temp_num, temp_node[0], temp_node[1], temp_node[2], node_line, K);
		}

		//将每个点属于哪些单元赋值
		for (ULL j = 0; j < 3; j++) {
			node_line[temp_node[j]].parentElement.insert(temp_num);
			//node_line[temp_node[j]].parentElement[node_line[temp_node[j]].parentElement_num] = temp_num;
			//node_line[temp_node[j]].parentElement_num++;
		}
	}

	//input_materials

	//材料属性数据规定E、miu、单元起点、单元终点
	double E_tmp, miu_tmp;
	ULL start, end;
	ULL material_num;
	inputfile >> material_num;
	for (ULL i = 1; i <= material_num; i++) {
		inputfile >> E_tmp >> miu_tmp >> start >> end;
		if (start == end && start != 0) {
			//单元属性
			element_line[start].set_E_miu(E_tmp, miu_tmp);
		}
		else if (start == 0) {
			for (ULL j = 1; j <= element_num; j++) {
				element_line[j].set_E_miu(E_tmp, miu_tmp);
			}
		}
		else {
			//单元组属性
			for (ULL j = start; j <= end; j++) {
				element_line[j].set_E_miu(E_tmp, miu_tmp);
			}
		}
	}
	inputfile.close();

	return true;
}
bool Test::draw_element() {
	QPolygonF points;
	QPen pen = QPen(Qt::black);
	pen.setColor(Qt::black);
	for (int i = 0; i < element_num; i++) {
		tri_item[i].setPen(pen);
		tri_item[i].setBrush(QColor(255, 255, 255));
		tri_item[i].setFlag(QGraphicsItem::ItemIsSelectable, true);
		for (int j= 0; j < 3; j++) 
			points.append(QPointF(element_line[i + 1].x[j] * EXPEND_COEF, -element_line[i + 1].y[j] * EXPEND_COEF));
		tri_item[i].setPolygon(points);
		tri_item[i].element_id = i + 1;
		scene->addItem(&tri_item[i]);
		points.clear();
	}
	return true;
}
bool Test::draw_node() {
	QPen pen = QPen(Qt::black);
	pen.setColor(Qt::black);
	pen.setWidth(1);
	for (int i = 1; i <= node_num; i++) {
		node_line[i].setPen(pen);
		node_line[i].setBrush(QColor(255, 255, 255, 80));
		node_line[i].setFlag(QGraphicsItem::ItemIsSelectable, true);
		node_line[i].setRect(QRectF((node_line[i].X * EXPEND_COEF - R*0.5), -(node_line[i].Y * EXPEND_COEF + R*0.5), R, R));
		scene_behind->addItem(&node_line[i]);
	}
	return true;
}
void Test::confirm_click() {
	QString tmp;
	ui.textBrowser->insertPlainText("计算中！\n");
	clock_t start, finish;
	double total;
	start = clock();
	if (if_input_model && if_input_boundary) {
		if (compute_K_Delta_F(element_num, node_num, constraint_num, constraint_level, K_resize, F_resize, Delta, element_line)) {
			ui.textBrowser->insertPlainText("计算成功！\n");
			redraw();
			ui.textBrowser->insertPlainText("重绘成功！\n");
		}
		else
			ui.textBrowser->insertPlainText("计算失败！\n");
	}
	finish = clock();
	total = (double)(finish - start) / CLOCKS_PER_SEC;
	ui.textBrowser->insertPlainText(tmp.asprintf("CPU占用的总时间：%fs\n", total));
	ui.groupBox->show();
}
void Test::clean_click() {
	scene->removeItem(tri_item);
	scene->update();
	//释放内存
	free(Delta);
	free(F);
	free(F_resize);
	free(K);
	free(K_resize);
	//释放内存
	delete[] node_line, element_line, constraint_line, constraint_level, constraint_line;
	delete[] tri_item, tri_item_new;//, cons_item, force_item;
}
void Test::close_click() {
	Test::close();
}
bool Test::redraw() {
	
	scene_new = new QGraphicsScene(this);
	tri_item_new = new GraphicsTriItem[element_num];
	
	double positve_sum[6] = { 0 };
	ULL positve_num[6] = { 0 };
	double negtive_sum[6] = { 0 };
	ULL negtive_num[6] = { 0 };
	for (ULL i = 1; i <= element_num; i++) {
		double* element_epsilon_tmp = element_line[i].epsilon;
		double* element_sigma_tmp = element_line[i].sigma;
		for (ULL j = 0; j < 3; j++) {
			if (element_epsilon_tmp[j] >= 0) {
				positve_sum[j] += element_epsilon_tmp[j];
				positve_num[j]++;
			}
			else {
				negtive_sum[j] += element_epsilon_tmp[j];
				negtive_num[j]++;
			}

			if (element_sigma_tmp[j] >= 0) {
				positve_sum[j + 3] += element_sigma_tmp[j];
				positve_num[j + 3]++;
			}
			else {
				negtive_sum[j + 3] += element_sigma_tmp[j];
				negtive_num[j + 3]++;
			}
		}
	}
	for (ULL i = 0; i < 6; i++) {
		positive_mean[i] = positve_sum[i] / positve_num[i];
		negtive_mean[i] = negtive_sum[i] / negtive_num[i];
	}
	for (ULL i = 1; i <= element_num; i++) {
		double* element_epsilon_tmp = element_line[i].epsilon;
		double* element_sigma_tmp = element_line[i].sigma;
		for (ULL j = 0; j < 3; j++) {
			if (element_epsilon_tmp[j] >= 0)
				positive_std[j] += pow(element_epsilon_tmp[j] - positive_mean[j], 2);
			else
				negtive_std[j] += pow(element_epsilon_tmp[j] - negtive_mean[j], 2);

			if (element_sigma_tmp[j] >= 0)
				positive_std[j + 3] += pow(element_sigma_tmp[j] - positive_mean[j + 3], 2);
			else
				negtive_std[j + 3] += pow(element_sigma_tmp[j] - negtive_mean[j + 3], 2);
		}
	}
	for (ULL i = 0; i < 6; i++) {
		positive_std[i] = sqrt(positive_std[i]);
		negtive_std[i] = sqrt(negtive_std[i]);
	}
	
	QPolygonF points;
	QPen pen = QPen(Qt::black);
	int e_s_offset = (is_e_s ? 0 : 3) + select_kind;
	for (ULL i = 0; i < element_num; i++) {
		double *measure_tmp = (is_e_s ? element_line[i + 1].epsilon : element_line[i + 1].sigma);
		double measure_k;// = measure_tmp[select_kind] * measure_len[p_n_offsets];
		if (measure_tmp[select_kind] > 0) 
			measure_k = (measure_tmp[select_kind] - positive_mean[e_s_offset]) / (constrast * concentrate_rate * positive_std[e_s_offset]);
		else
			measure_k = -(measure_tmp[select_kind] - negtive_mean[e_s_offset]) / (constrast * concentrate_rate * negtive_std[e_s_offset]);
		if (fabs(measure_k) < 1)
			measure_k += 1;////////
		else
			measure_k = 2;
		
		QColor color;
		measure_k *= 127;
		color.setHsv(((measure_tmp[select_kind] > 0) ? 0 : 240), ((fabs(measure_k)<0.01)?0:ceil(measure_k)), toumingdu);
		pen.setColor((boundary_color) ? Qt::black : color);
		tri_item_new[i].setPen(pen);
		tri_item_new[i].setBrush(color);
		tri_item_new[i].setFlag(QGraphicsItem::ItemIsSelectable, true);
		for (int j = 0; j < 3; j++)
			points.append(QPointF(element_line[i + 1].x[j] * EXPEND_COEF, -element_line[i + 1].y[j] * EXPEND_COEF));
		tri_item_new[i].setPolygon(points);
		tri_item_new[i].element_id = i + 1;
		scene_new->addItem(&tri_item_new[i]);
		points.clear();
	}
	
	ui.graphicsView->setScene(scene_new);
	return true;
}
bool Test::select_redraw() {
	QPen pen = QPen(Qt::black);
	int e_s_offset = (is_e_s ? 0 : 3) + select_kind;
	for (ULL i = 0; i < element_num; i++) {
		/*double* measure_tmp = (is_e_s ? element_line[i + 1].epsilon : element_line[i + 1].sigma);
		int p_n_offsets = ((measure_tmp[select_kind] > 0) ? 0 : 3) + e_s_offset;
		double measure_k = measure_tmp[select_kind] * measure_len[p_n_offsets];

		QColor color;*/
		double* measure_tmp = (is_e_s ? element_line[i + 1].epsilon : element_line[i + 1].sigma);
		double measure_k;// = measure_tmp[select_kind] * measure_len[p_n_offsets];
		if (measure_tmp[select_kind] > 0)
			measure_k = (measure_tmp[select_kind] - positive_mean[e_s_offset]) / (constrast * concentrate_rate * positive_std[e_s_offset]);
		else
			measure_k = -(measure_tmp[select_kind] - negtive_mean[e_s_offset]) / (constrast * concentrate_rate * negtive_std[e_s_offset]);
		if (fabs(measure_k) < 1)
			measure_k += 1;////////
		else
			measure_k = 2;
		


		QColor color;
		measure_k *= 127;
		color.setHsv(((measure_tmp[select_kind] > 0) ? 0 : 240), ((fabs(measure_k) < 0.01) ? 0 : ceil(measure_k)), toumingdu);
		pen.setColor((boundary_color) ? Qt::black : color);
		tri_item_new[i].setPen(pen);
		tri_item_new[i].setBrush(color);
	}
	return true;
}
void Test::color_change() {
	constrast = ui.horizontalSlider->value();
	toumingdu = ui.horizontalSlider_2->value();
	boundary_color = ui.checkBox->isChecked();
	int tmp_num = block1->checkedId();
	if (tmp_num >= 3) {
		is_e_s = true;
		select_kind = tmp_num - 3;
	}
	else {
		is_e_s = false;
		select_kind = tmp_num;
	}
	select_redraw();
	
}
void Test::get_two_number() {
	double check_x = ui.lineEdit->text().toDouble();
	double check_y = ui.lineEdit_3->text().toDouble();
	QPoint check_point(EXPEND_COEF * check_x, -EXPEND_COEF * check_y);

	QList<QGraphicsItem*> element_items;//碰撞到的单元
	QList<QGraphicsItem*> node_items;//碰撞到的节点
	element_items = scene->items(check_point);
	node_items = scene_behind->items(check_point);
	check_point_information(check_x, check_y,element_items, node_items);
}
void Test::check_point_information(double x,double y,QList<QGraphicsItem*> element_items, QList<QGraphicsItem*> node_items) {
	double check_x = x;
	double check_y = y;
	QPoint check_point(EXPEND_COEF * check_x, -EXPEND_COEF * check_y);
	element_items = scene->items(check_point);
	node_items = scene_behind->items(check_point);
	QString tmp;
	
	//寻找最近节点距离
	double min_distance = 1000000;
	int min_node_index = 0;
	double tmp_distance = 0;
	ULL select_node_id = 0;
	ULL select_element_id = 0;
	
	
	if (!node_items.isEmpty()) {
		Node* select_node_item = (Node*)node_items[0];
		select_node_id = select_node_item->nodenum;
		ui.textBrowser->insertPlainText(tmp.asprintf("靠近该点处有一节点(%.2f,%.2f)\n", select_node_item->X, select_node_item->Y));
		ui.textBrowser->insertPlainText(tmp.asprintf("其节点位移为(%.2f,%.2f),其父节点有：\n", Delta[2 * (select_node_id - 1)], Delta[2 * (select_node_id - 1) + 1]));
		std::set<ULL>::iterator i;
		for (i = select_node_item->parentElement.begin(); i != select_node_item->parentElement.end(); i++) {
			ui.textBrowser->insertPlainText(tmp.asprintf("%d\n", *i));
		}
	}
	if (!element_items.isEmpty()) {
		GraphicsTriItem* select_element_item = (GraphicsTriItem*)element_items[0];
		select_element_id = select_element_item->element_id;
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在一单元%d内\n", select_element_id));
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在单元内位移为(%.2f,%.2f)\n", element_line[select_element_id].get_u(check_x, check_y), element_line[select_element_id].get_v(check_x, check_y)));
		ui.textBrowser->insertPlainText(tmp.asprintf("该单元内应变\n"));
		ui.textBrowser->insertPlainText(tmp.asprintf("ε_x:%.4f\n", element_line[select_element_id].epsilon[0]));
		ui.textBrowser->insertPlainText(tmp.asprintf("ε_y:%.4f\n", element_line[select_element_id].epsilon[1]));
		ui.textBrowser->insertPlainText(tmp.asprintf("γ_xy:%.2f\n", element_line[select_element_id].epsilon[2]));
		ui.textBrowser->insertPlainText(tmp.asprintf("该单元内应力\n"));
		ui.textBrowser->insertPlainText(tmp.asprintf("σ_x:%.2f\n", element_line[select_element_id].sigma[0]));
		ui.textBrowser->insertPlainText(tmp.asprintf("σ_y:%.2f\n", element_line[select_element_id].sigma[1]));
		ui.textBrowser->insertPlainText(tmp.asprintf("τ_xy:%.2f\n", element_line[select_element_id].sigma[2]));
	}
	if(node_items.isEmpty()&& element_items.isEmpty())
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在模型外\n"));
	ui.textBrowser->insertPlainText("\n");
}
void Test::get_QPointF(QPointF point) {
	double check_x = point.x() / EXPEND_COEF;
	double check_y = -point.y() / EXPEND_COEF;
	QPoint check_point(EXPEND_COEF * check_x, -EXPEND_COEF * check_y);
	
	QString tmp;
	ui.textBrowser->insertPlainText(tmp.asprintf("双击点:(%.2f,%.2f)\n", check_x, check_y));
	
	QList<QGraphicsItem*> element_items;//碰撞到的单元
	element_items = scene->items(point);
	check_QPointF_information(check_x,check_y,element_items);
}
void Test::check_QPointF_information(double x,double y,QList<QGraphicsItem*> element_items) {
	QString tmp;
	double check_x = x;
	double check_y = y;
	if (!element_items.isEmpty()) {
		GraphicsTriItem* select_element_item = (GraphicsTriItem*)element_items[0];
		ULL select_element_id = select_element_item->element_id;
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在一单元%d内\n", select_element_id));
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在单元内位移为(%.2f,%.2f)\n", element_line[select_element_id].get_u(check_x, check_y), element_line[select_element_id].get_v(check_x, check_y)));
		ui.textBrowser->insertPlainText(tmp.asprintf("该单元内应变\n"));
		ui.textBrowser->insertPlainText(tmp.asprintf("ε_x:%.4f\n", element_line[select_element_id].epsilon[0]));
		ui.textBrowser->insertPlainText(tmp.asprintf("ε_y:%.4f\n", element_line[select_element_id].epsilon[1]));
		ui.textBrowser->insertPlainText(tmp.asprintf("γ_xy:%.4f\n", element_line[select_element_id].epsilon[2]));
		ui.textBrowser->insertPlainText(tmp.asprintf("该单元内应力\n"));
		ui.textBrowser->insertPlainText(tmp.asprintf("σ_x:%.2f\n", element_line[select_element_id].sigma[0]));
		ui.textBrowser->insertPlainText(tmp.asprintf("σ_y:%.2f\n", element_line[select_element_id].sigma[1]));
		ui.textBrowser->insertPlainText(tmp.asprintf("τ_xy:%.2f\n", element_line[select_element_id].sigma[2]));
	}else
		ui.textBrowser->insertPlainText(tmp.asprintf("该点在模型外\n"));
	ui.textBrowser->insertPlainText("\n");
}
Test::~Test() {
	delete[] tri_item, tri_item_new; //,cons_item, force_item;
	delete[] node_line, element_line, constraint_line, constraint_level;
	delete scene,scene_behind,scene_new, block1;
	free(F);
	free(K);
	free(Delta);
	free(F_resize);
	free(K_resize);
}