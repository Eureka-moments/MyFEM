#include "compute.h"
/*预定义的输入数据
* ！！！！规定！！！！！
*1.全部计算数据使用double类型双精度
*2.按需求增加中间计算留存数据（可变）
*3.规范统一变量名称，方便调试
*！！！！！！！！！！！！*
*/

//任务目标实现
//1.计算外力分布向量
// （坐标与节点/单元的对应问题。）
//2.计算单元刚度矩阵
//3.合并计算总刚矩阵
//4.求解总刚的位移-外力大矩阵线性方程
//5.求得位移后回代得到其他一些力学分析信息

bool update_original_delta(ULL element_num,ULL constraint_num, ULL* constraint_level,double* F_resize, double* Delta, Element* element_line) {
    /*for (ULL i = 0; i < constraint_num; i++) {
        Delta[constraint_line[i]] = 0;
    }*/
    double* delta_resize_tmp = F_resize;
    double* delta_tmp = Delta;
    ULL cpy_num = 0;
    //根据算出来的delta_resize更新原Delta
    for (ULL i = 1; i <= constraint_num + 1; i++) {//根据标记点偏移量进行memcpy
        cpy_num = constraint_level[i];
        if (cpy_num != 0)
            memcpy(delta_tmp, delta_resize_tmp, cpy_num * SIZE_DOUBLE);
        delta_resize_tmp += cpy_num;
        delta_tmp += (cpy_num + 1);
    }

    for (ULL i = 1; i <= element_num; i++) {
        element_line[i].set_uv(Delta);
    }
    return true;
}
bool update_resize_F(ULL constraint_num, ULL* constraint_level,double* F,double* F_resize) {
    double* F_resize_tmp = F_resize;
    double* F_tmp = F;
    ULL cpy_num = 0;
    
    //创建resize_F
    for (ULL i = 1; i <= constraint_num + 1; i++)
    { //根据标记点偏移量进行memcpy
        cpy_num = constraint_level[i];
        if (cpy_num != 0)
            memcpy(F_resize_tmp, F_tmp, cpy_num * sizeof(double));
        F_resize_tmp += cpy_num;
        F_tmp += cpy_num + 1;
    }

    return true;
}
bool update_resize_K(ULL constraint_num, ULL* constraint_line, ULL* constraint_level, double* K,double* K_resize) {
    double* K_resize_tmp = K_resize; //存储指针
    double* K_tmp = K;               //被cpy指针
    ULL cpy_num = 0;  // cpy长度

    ULL level_seg_num; //每层内段数
    ULL advanced_num = 1;    //层内自增数代表层内行数增加
    ULL belong_cons;   //当前段属于哪一约束
    for (ULL i = 1; i <= constraint_num + 1; i++)
    { //循环每层
        if (constraint_level[i] != 0)
        { //若去除行两行相邻则只需要偏移而不需要cpy
            if (i == 1)
                level_seg_num = 1;
            else
                level_seg_num = constraint_level[i] * i;
            for (ULL j = 1; j <= level_seg_num; j++)
            {                        //循环每段
                belong_cons = j % i; //当前段属于哪一约束
                if (i == 1)          //如果是第一层则初始重置为三角
                    cpy_num = (constraint_level[i] + 1) * constraint_level[i] / 2;
                else
                {                         //不是第一层
                    if (belong_cons == 0) //是每行最后一段
                    {
                        cpy_num += advanced_num; //叠加行尾偏移
                        advanced_num++;
                        if (j != level_seg_num) //如果已经是该层最后一段
                            continue;           //还有下一行则叠加到下一行第一段cpy
                    }                           //每段长度，即cpy长度cpy_num
                    else                        //每层除尾外叠加固定长度
                        cpy_num += constraint_level[belong_cons];
                }
                if (cpy_num != 0)
                    memcpy(K_resize_tmp, K_tmp, cpy_num * sizeof(double));
                K_resize_tmp += cpy_num; //指针偏移增加量
                K_tmp += cpy_num + 1;    //指针偏移增加量（+1由于标记/洞占位）
                cpy_num = 0;
            }
            K_tmp -= 1;
        }
        advanced_num = 1;                //每层结束将自增数重置
        K_tmp += constraint_line[i] + 1; //每层结束标记行占位跳过
    }
    return true;
}
bool compute_K_Delta_F(ULL element_num, ULL node_num, ULL constraint_num, ULL* constraint_level, double* K_resize, double* F_resize, double* Delta, Element* element_line) {
    int *ipiv = (int*)calloc(2 * node_num - constraint_num, sizeof(int));
	
    LAPACKE_dsptrf(LAPACK_ROW_MAJOR, 'L', 2 * node_num - constraint_num, K_resize, ipiv);
    LAPACKE_dsptrs(LAPACK_ROW_MAJOR, 'L', 2 * node_num - constraint_num, 1, K_resize, ipiv, F_resize, 1);
	
    free(ipiv);
    return update_original_delta(element_num, constraint_num, constraint_level, F_resize, Delta, element_line);
}
double vector_dot_product(double* a, double* b, int n) { return cblas_ddot(n, a, 1, b, 1); }
void add_vector(double* a, double* b, double* c, int n) {
    cblas_daxpy(n, 1, b, 1, c, 1);
    cblas_daxpy(n, 1, a, 1, c, 1);
}
void subtract_vector(double* a, double* b, double* c, int n) {
    cblas_daxpy(n, -1, b, 1, c, 1);
    cblas_daxpy(n, 1, a, 1, c, 1);
}
void multiply_scalar(double s, double* a, int n) { cblas_dscal(n, s, a, 1); }
double vecor_length(double* a, double* b, int n)
{
    double* c = (double*)calloc(n, sizeof(double));
    subtract_vector(a, b, c, n);
    double result = cblas_dnrm2(n, c, 1);
    free(c);
    return result;
}
