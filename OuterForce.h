#pragma once

//定义外力包括集中力和边缘分布力
class Force {
public:
    Force() {};
    ~Force() {};
    //作用方向向量
    double x;
    double y;
    //作用点
    double x_s;
    double y_s;
    double x_e;
    double y_e;
    //作用大小
    double f;


private:

};

//定义外力体力给出了加速度大小和方向
class BodyForce {
public:
    BodyForce(double x, double y, double a) : x(x), y(y), a(a) {};
    ~BodyForce() {};
    //作用方向向量
    double x;
    double y;
    //作用大小
    double a;
private:

};
