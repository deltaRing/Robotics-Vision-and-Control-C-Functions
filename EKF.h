#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

// 原作者：
//https://blog.csdn.net/zhangquan2015/article/details/79264540

using namespace std;
using namespace Eigen;

MatrixXd EKF1() {
	double T = 1; // 雷达扫描周期
	double N = 60.0 / T; // 总的采样次数
	MatrixXd X = MatrixXd::Zero(4, N); // 目标初始位置、速度
	X.col(0) << -100, 2, 200, 20; // 
	MatrixXd Z = MatrixXd::Zero(1, N); // 传感器对位置的观测

	double delta_w = 1e-3; // 真实航迹的调整值
	MatrixXd Q(2, 2);
	Q << 0.5, 0, 0, 1;
	Q = Q * delta_w; // 过程噪声方差
	MatrixXd G(4, 2);
	G << pow(T, 2) / 2, 0, T, 0, 0, pow(T, 2) / 2, 0, T; // 过程噪声驱动矩阵
	double R = 5; // 观测噪声方差
	MatrixXd F(4, 4);
	F << 1, T, 0, 0, 0, 1, 0, 0, 0, 0, 1, T, 0, 0, 0, 1; // 状态转移矩阵
	double x0 = 200; // 观察位置
	double y0 = 300;
	MatrixXd Xstation(2, 1);
	Xstation << x0, y0;

	for (int t = 1; t < N; t++) {
		X.col(t) = F * X.col(t - 1) + G * sqrt(Q.array()).matrix() * MatrixXd::Random(2, 1); // 目标真实轨迹
	}

	for (int t = 0; t < N; t++) {
		Z(t) = sqrt(pow(X.col(t)(0) - Xstation(0), 2) + pow(X.col(t)(1) - Xstation(1), 2)) + sqrt(R) * 0.1; // 原本应该是randn
	}

	// ----------------------------------------------- EKF ----------------------------------------------------------------
	MatrixXd Xekf = MatrixXd::Zero(4, N); // 目标初始位置、速度
	Xekf.col(0) = X.col(0); //卡尔曼滤波状态初始化
	MatrixXd P0 = MatrixXd::Identity(4, 4); // 协方差初始化

	for (int i = 1; i < N; i++) {
		MatrixXd Xn = F * Xekf.col(i - 1); // 预测
		MatrixXd P1 = F * P0 * F.transpose() + G * Q * G.transpose(); // 预测误差协方差
		double dd = sqrt(pow(X.col(i)(0) - Xstation(0), 2) + pow(X.col(i)(1) - Xstation(1), 2)); // 观测预测
		// 求雅可比矩阵H
		MatrixXd H(1, 4); // 求一阶近似
		H << (Xn(0, 0) - x0) / dd, 0, (Xn(2, 0) - y0) / dd, 0;
		MatrixXd K = P1 * H.transpose() * ((H * P1 * H.transpose()).array() + R).inverse().matrix(); // 增益
		Xekf.col(i) = Xn + K * (Z.col(i).array() - dd).matrix();
		P0 = (MatrixXd::Identity(4, 4) - K * H) * P1; // 滤波误差协方差更新
	}
	return Xekf;
}
