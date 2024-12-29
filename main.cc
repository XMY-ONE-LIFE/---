#include <iostream>                                  // 包含输入输出流库
#include<Eigen/Dense>                               // 包含Eigen库的Dense模块,用于矩阵运算
#include<cmath>                                     // 包含数学函数库
#include <tuple>                                    // 包含元组库
#include<algorithm>                                 // 包含算法库
#include<vector>                                    // 包含向量容器库
#include <thread>                                   // 包含线程库
#include <matplot/matplot.h>                        // 包含Matplot++绘图库
#include <set>                                      // 包含集合容器库

using namespace std;                                // 使用标准命名空间
using namespace Eigen;                              // 使用Eigen命名空间
struct PathPoint {                                  // 定义路径点结构体
    double x;                                       // x坐标
    double y;                                       // y坐标
};

void draw_single_line(PathPoint &p1,PathPoint &p2)
{
    std::vector<double> x, y;                       // 定义存储路径点的向量
    
    Eigen::VectorXd X(6);                    // 定义6维向量存储x方向的约束条件
    Eigen::VectorXd Y(6);                    // 定义6维向量存储y方向的约束条件
 
    double x0 = p1.x;                      // 起点x坐标
    double dx0 = 5.0;                              // 起点x方向速度
    double ddx0 = 0.0;                             // 起点x方向加速度
    double x1 = p2.x;                       // 终点x坐标
    double dx1 = 5.0;                              // 终点x方向速度
    double ddx1 = 0.0;                             // 终点x方向加速度
 
    double y0 = p1.y;                      // 起点y坐标
    double dy0 = 0.0;                              // 起点y方向速度
    double ddy0 = 0.0;                             // 起点y方向加速度
    double y1 = p2.y;                       // 终点y坐标
    double dy1 = 0.0;                              // 终点y方向速度
    double ddy1 = 0.0;                             // 终点y方向加速度
 
    X << x0, dx0, ddx0, x1, dx1, ddx1;            // 设置x方向约束条件向量
    Y << y0, dy0, ddy0, y1, dy1, ddy1;            // 设置y方向约束条件向量
 
    double t0 = 0;                                 // 起始时间
    double t1 = 3;                                 // 终止时间
 
    Eigen::MatrixXd T(6, 6);                       // 定义6x6矩阵存储时间多项式系数
 
 
    T << 1, t0, pow(t0, 2), pow(t0, 3), pow(t0, 4), pow(t0, 5),           // 设置时间多项式系数矩阵第1行
        0, 1, 2 * t0, 3 * pow(t0, 2), 4 * pow(t0, 3), 5 * pow(t0, 4),     // 设置时间多项式系数矩阵第2行
        0, 0, 2, 6 * t0,12 * pow(t0, 2), 20 * pow(t0, 3),                 // 设置时间多项式系数矩阵第3行
        1, t1, pow(t1, 2), pow(t1, 3), pow(t1, 4), pow(t1, 5),            // 设置时间多项式系数矩阵第4行
        0, 1, 2 * t1, 3 * pow(t1, 2), 4 * pow(t1, 3),5 * pow(t1, 4),      // 设置时间多项式系数矩阵第5行
        0, 0, 2, 6 * t1, 12 * pow(t1, 2), 20 * pow(t1, 3);                // 设置时间多项式系数矩阵第6行
 
    //计算A和B两个系数矩阵
    Eigen::MatrixXd A = T.inverse() * X;           // 计算x方向多项式系数
    Eigen::MatrixXd B = T.inverse() * Y;           // 计算y方向多项式系数
    
    //保存横纵向坐标和横纵向速度用于画图
    vector<double>x_,y_,v_lon,v_lat;               // 定义存储轨迹点和速度的向量
 
    vector<double>time;                            // 定义存储时间点的向量
    int count = 0;                                 // 计数器初始化
    for(double t=t0;t<t1+0.05;t+=0.05) {          // 以0.05s为间隔生成时间序列
        count++;                                   // 计数器增加
        time.push_back(t);                         // 将时间点加入向量
 
    }
 
    MatrixXd temp1(1,6),temp2(1,6);               // 定义临时矩阵(1行6列)用于计算位置和速度 
    for(int i=0;i<count;i++){                     // 遍历所有时间点
 
        temp1<<1,time[i],pow(time[i], 2),pow(time[i], 3),pow(time[i], 4),pow(time[i], 5);    // 设置位置计算的时间多项式系数
        x_.push_back((temp1*A)(0,0));             // 计算并存储x坐标
        y_.push_back((temp1*B)(0,0));             // 计算并存储y坐标
 
        temp2<<0,1,2*time[i],3*pow(time[i],2),4*pow(time[i],3),5*pow(time[i],4);             // 设置速度计算的时间多项式系数
        v_lon.push_back((temp2*A)(0,0));          // 计算并存储纵向速度
        v_lat.push_back((temp2*B)(0,0));          // 计算并存储横向速度
    }

     
    matplot::hold(matplot::on);// 设置绘图保持（hold）状态，这样后续绘制的曲线会添加到已有的图形上，而不会覆盖之前的
    
    matplot::plot(x_, y_);                        // 绘制轨迹图
 
   
    
    // matplot::plot(time,v_lon);                 // 绘制纵向速度图
    // matplot::plot(time,v_lat);                 // 绘制横向速度图
 
  
}


 
int main () {                                       // 主函数入口

    std::vector<PathPoint> startPoints = {// 创建多个起点和终点的示例，这里简单创建了3组点对，你可以根据实际需求增加或修改
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75},
        {0.0, -1.75}
    };
    std::vector<PathPoint> endPoints = {
        {20.0, 1.75},
        {20.0, 1.35},
        {20.0, 0.95},
        {20.0, 0.55},
        {20.0, 0.15},
        {20.0, -0.35},
        {20.0, -0.75},
        {20.0, -1.15},
        {20.0, -1.55},
        {20.0, -1.95},
        {20.0, -2.35},
        {20.0, -2.75}
    };

    
    if (startPoints.size()!= endPoints.size()) {// 确保起点和终点的数量是一致的，否则逻辑会出错
        std::cerr << "起点和终点数量不一致，请检查输入数据！" << std::endl;
        return -1;
    }

   
    matplot::figure(1); // 创建图形窗口1，用于绘制多条轨迹曲线

    
    for (size_t i = 0; i < startPoints.size(); ++i) {// 循环遍历每一组起点和终点，调用绘制函数绘制曲线
        draw_single_line(startPoints[i], endPoints[i]);
    }

    
    matplot::show();// 显示所有图形（所有曲线都会展示在图形窗口1中）

    return 0;

}