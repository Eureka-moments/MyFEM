#include "element.h"

ULL Element::get_line_index(ULL i, ULL j)
{
    if (i >= j)
        return i * (i - 1) / 2 + j - 1;
    else
        return j * (j - 1) / 2 + i - 1;;
}

//单元信息处理
void Element::set(ULL num, ULL node1, ULL node2, ULL node3, Node* nodeline_input,double* K_input) {
    elementnum = num;
    node[0] = node1;
    node[1] = node2;
    node[2] = node3;
    //K矩阵改变
    K = K_input;
    //node_line
    node_line = nodeline_input;
    //从读入节点信息初始化单元内节点坐标
    for (int i = 0; i < NODE_NUM_IN_ELEMENT; i++) {
        x[i] = node_line[node[i]].X;
        y[i] = node_line[node[i]].Y;
    }
    //xy改变
    is_xy = true;
    //更新
    check_update();
}

void Element::set_uv(double* Delta) {//节点位移改变(产生)后修改对应单元的内部数值
    for (int i = 0; i < NODE_NUM_IN_ELEMENT; i++) {
        u[i] = Delta[2 * (node[i] - 1)];
        v[i] = Delta[2 * (node[i] - 1) + 1];
        //更新xy坐标
        x[i] += u[i];
        y[i] += v[i];
        //更新xy序列
        node_line[node[i]].X = x[i];
        node_line[node[i]].Y = y[i];
    }
    //xy改变
    is_xy = true;
    //uv改变
    is_uv = true;

    //更新
    check_update();
}
void Element::set_E_miu(double E_tmp, double miu_tmp) {
    E = E_tmp;
    miu = miu_tmp;
    is_E_miu = true;
    check_update();
}
//获取单元内某点位移
double Element::get_u(double x_in, double y_in) {
    return alpha[0] + alpha[1] * x_in + alpha[2] * y_in;
}
double Element::get_v(double x_in, double y_in) {
    return beta[0] + beta[1] * x_in + beta[2] * y_in;
}
//private function
void Element::update_abcS() {
    int i, j, k;
    for (i = 0; i < NODE_NUM_IN_ELEMENT; i++) {
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        a[i] = x[j] * y[k] - x[k] * y[j];
        b[i] = y[j] - y[k];
        c[i] = x[k] - x[j];
    }
    S = 0.5 * (y[2] * (x[1] - x[0]) - y[1] * (x[2] - x[0]) - y[0] * (x[1] - x[2]));
    S_scala = 0.5 * (1 / S);
}
void Element::update_epsilon() {
    //位移偏导矩阵——应变量
    //double epsilon_temp[PROBLEM_DIMENSION * PROBLEM_DIMENSION] = { 0 };//4
    
    epsilon[0] = alpha[1];
    epsilon[1] = alpha[5];
    epsilon[2] = alpha[2] + alpha[4];
    /*cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2, 2, 3, S_scala, b, 3, u, 2, 0, epsilon_temp, 2);
    epsilon[0] = epsilon_temp[0];
    epsilon[1] = epsilon_temp[3];
    epsilon[2] = epsilon_temp[1] + epsilon_temp[2];*/
}

void Element::update_sigma() {
    sigma[0] = E_scala * (epsilon[0] + miu * epsilon[1]);
    sigma[1] = E_scala * (epsilon[1] + miu * epsilon[0]);
    sigma[2] = 0.5 * E_scala * (1 - miu) * epsilon[2];
}

void Element::update_E_scala() {
    E_scala = E / (1 - miu * miu);
}

void Element::update_alpha() {
    //更新alpha用来计算应力应变
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2, 3, 3, S_scala, u, 3, a, 3, 0, alpha, 3);
}
void Element::update_K() {

    double temp_k = E * t / (4*(1-miu*miu)* S);
    double miu_tmp = temp_k * (1 - miu)/2;

    double matrix[2][2] = { 0 };
    double bb_tmp = 0;
    double cc_tmp = 0;
    double bc_tmp = 0;
    double cb_tmp = 0;

    for (ULL i = 0; i < 3; i++) {
        for (ULL j = 0; j <= i; j++) {
            bb_tmp = b[i] * b[j];
            cc_tmp = c[i] * c[j];
            bc_tmp = b[i] * c[j];
            cb_tmp = c[i] * b[j];
            matrix[0][0] = temp_k * bb_tmp + miu_tmp * cc_tmp;
            matrix[0][1] = temp_k * miu * bc_tmp + miu_tmp * cb_tmp;
            matrix[1][1] = temp_k * cc_tmp + miu_tmp * bb_tmp;
			matrix[1][0] = temp_k * miu * cb_tmp + miu_tmp * bc_tmp;
            //先给主对角线赋值
            if (i == j) {
                for (ULL m = 1; m <= 2; m++)
                    for (ULL n = 1; n <= m; n++)
                        K[get_line_index(2 * (node[i]-1) + m, 2 * (node[j]-1) + n)] += matrix[m - 1][n - 1];
            }
            else {//然后是其他部分
                for (ULL m = 1; m <= 2; m++)
                    for (ULL n = 1; n <= 2; n++)
                        K[get_line_index(2 * (node[i]-1) + m, 2 * (node[j]-1) + n)] += matrix[m - 1][n - 1];
            }
        }
    }
}


bool Element::check_update() {
    if (is_uv || is_xy) {//对于xy和uv的更新有一些部分是要一起更新的，而只有E、miu更新则情况不同
        if (is_xy) {
            update_abcS();
            is_xy = false;
        }
        update_alpha();
        update_epsilon();
        update_E_scala();
        update_sigma();
        if (is_E_miu) {
            update_K();
            is_E_miu = false;
        }
        is_uv = false;
        return true;
    }
    if (is_E_miu) {
        update_E_scala();
        update_sigma();
        update_K();
        is_E_miu = false;
        return true;
    }
    return false;//返回false是表示没有任何更新
}
