//#include "inputfile.h"
//
//
//extern ULL element_num;
//extern ULL node_num;//节点总数(考虑写成全局)
//extern ULL constraint_num;//约束总数(考虑写成全局)
//extern Node* node_line;
//extern Element* element_line;
//extern ULL* constraint_line;
//extern double* Delta;
//extern double* Delta_resize;
//extern double* F;
//extern double* F_resize;
//extern double* K;
//extern double* K_resize;
////extern double zou;
////extern double t;
//
//bool Input_Node(char* node_location) {
//
//    ULL temp_num;
//    double temp_x, temp_y;
//    std::ifstream inputfile;
//    inputfile.open(node_location);
//    inputfile >> node_num;
//    node_line = new Node[node_num + 1];
//    Delta = (double*)mkl_calloc(2 * node_num, SIZE_DOUBLE, 64);
//    memset(Delta, 0, 2 * node_num * SIZE_DOUBLE);
//    F = (double*)mkl_calloc(2 * node_num, SIZE_DOUBLE, 64);
//    K = (double*)mkl_calloc(((2 * node_num + 1) * node_num), SIZE_DOUBLE, 64);
//    memset(K, 0, ((2 * node_num + 1) * node_num) * SIZE_DOUBLE);
//    //初始化点
//    for (ULL i = 1; i <= node_num; i++) {
//        inputfile >> temp_num >> temp_x >> temp_y;
//        node_line[i].set(temp_num, temp_x, temp_y);
//    }
//    inputfile.close();
//    return true;
//}
//bool Input_Element(char* element_location) {
//    ULL temp_node[3];
//    ULL temp_num;
//    std::ifstream inputfile;
//    inputfile.open(element_location);
//
//    inputfile >> element_num;
//    element_line = new Element[element_num + 1];
//
//    //初始化单元
//    for (ULL i = 1; i <= element_num; i++) {
//        inputfile >> temp_num >> temp_node[0] >> temp_node[1] >> temp_node[2];
//        double x_1 = node_line[temp_node[1]].X - node_line[temp_node[0]].X;
//        double x_2 = node_line[temp_node[2]].X - node_line[temp_node[0]].X;
//        double y_1 = node_line[temp_node[1]].Y - node_line[temp_node[0]].Y;
//        double y_2 = node_line[temp_node[2]].Y - node_line[temp_node[0]].Y;
//        if (x_1 * y_2 - x_2 * y_1 < 0) {
//            element_line[i].set(temp_num, temp_node[0], temp_node[2], temp_node[1]);
//        }else {
//            element_line[i].set(temp_num, temp_node[0], temp_node[1], temp_node[2]);
//        }
//        
//        //将每个点属于哪些单元赋值
//        for (ULL j = 0; j < 3; j++) {
//            node_line[temp_node[j]].parentElement.insert(temp_num);
//            //node_line[temp_node[j]].parentElement[node_line[temp_node[j]].parentElement_num] = temp_num;
//            //node_line[temp_node[j]].parentElement_num++;
//        }
//    }
//    inputfile.close();
//
//    return true;
//}
//bool Input_Material(char* Material_location) {
//	
//
//	std::ifstream inputfile;
//	inputfile.open(Material_location);
//	//材料属性数据规定E、miu、单元起点、单元终点
//    double E_tmp, miu_tmp;
//    ULL start, end;
//	ULL material_num;
//	inputfile >> material_num;
//    for (ULL i = 1; i <= material_num; i++) {
//        inputfile >> E_tmp >> miu_tmp >> start >> end;
//        if (start == end && start != 0) {
//        //单元属性
//            element_line[start].set_E_miu(E_tmp, miu_tmp);
//		}
//        else if (start == 0) {
//            for (ULL j = 1; j <= element_num; j++) {
//                element_line[j].set_E_miu(E_tmp, miu_tmp);
//            }
//        }
//		else {
//			//单元组属性
//			for (ULL j = start; j <= end; j++) {
//                element_line[j].set_E_miu(E_tmp, miu_tmp);
//			}
//		}
//    }
//	inputfile.close();
//
//    return true;
//}
//bool Input_Model(char* model_location) {
//    //input_node
//    ULL temp_num;
//    double temp_x, temp_y;
//    std::ifstream inputfile;
//    inputfile.open(model_location);
//    inputfile >> node_num;
//	
//    node_line = new Node[node_num + 1];
//    Delta = (double*)calloc(2 * node_num, SIZE_DOUBLE);
//    F = (double*)calloc(2 * node_num, SIZE_DOUBLE);
//    K = (double*)calloc(((2 * node_num + 1) * node_num), SIZE_DOUBLE);
//    
//    //初始化点
//    for (ULL i = 1; i <= node_num; i++) {
//        inputfile >> temp_num >> temp_x >> temp_y;
//        node_line[i].set(temp_num, temp_x, temp_y);
//    }
//
//	
//
//	//input_element
//    ULL temp_node[3];
//    /*ULL temp_num;*/
//
//    inputfile >> element_num;
//    element_line = new Element[element_num + 1];
//
//    //初始化单元
//    for (ULL i = 1; i <= element_num; i++) {
//        inputfile >> temp_num >> temp_node[0] >> temp_node[1] >> temp_node[2];
//        double x_1 = node_line[temp_node[1]].X - node_line[temp_node[0]].X;
//        double x_2 = node_line[temp_node[2]].X - node_line[temp_node[0]].X;
//        double y_1 = node_line[temp_node[1]].Y - node_line[temp_node[0]].Y;
//        double y_2 = node_line[temp_node[2]].Y - node_line[temp_node[0]].Y;
//        if (x_1 * y_2 - x_2 * y_1 < 0) {
//            element_line[i].set(temp_num, temp_node[0], temp_node[2], temp_node[1]);
//        }
//        else {
//            element_line[i].set(temp_num, temp_node[0], temp_node[1], temp_node[2]);
//        }
//
//        //将每个点属于哪些单元赋值
//        for (ULL j = 0; j < 3; j++) {
//            node_line[temp_node[j]].parentElement.insert(temp_num);
//            //node_line[temp_node[j]].parentElement[node_line[temp_node[j]].parentElement_num] = temp_num;
//            //node_line[temp_node[j]].parentElement_num++;
//        }
//    }
//
//    //input_materials
//
//    //材料属性数据规定E、miu、单元起点、单元终点
//    double E_tmp, miu_tmp;
//    ULL start, end;
//    ULL material_num;
//    inputfile >> material_num;
//    for (ULL i = 1; i <= material_num; i++) {
//        inputfile >> E_tmp >> miu_tmp >> start >> end;
//        if (start == end && start != 0) {
//            //单元属性
//            element_line[start].set_E_miu(E_tmp, miu_tmp);
//        }
//        else if (start == 0) {
//            for (ULL j = 1; j <= element_num; j++) {
//                element_line[j].set_E_miu(E_tmp, miu_tmp);
//            }
//        }
//        else {
//            //单元组属性
//            for (ULL j = start; j <= end; j++) {
//                element_line[j].set_E_miu(E_tmp, miu_tmp);
//            }
//        }
//    }
//    inputfile.close();
//
//    return true;
//}