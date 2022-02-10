#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

// 房间温度在25摄氏度左右，测量误差为正负0.5摄氏度，方差0.25，R = 0.25。Q = 0.01, A = 1, T = 1, H = 1。
// 假定快时刻的温度值、测量值为23.9摄氏度，房间真实温度为24摄氏度，
//温度计在该时刻测量值为24.5摄氏度，偏差为0.4摄氏度。利用k - 1时刻温度值测量第k时刻的温度，其预计偏差为:
// P(k | k - 1) = P(k - 1) + Q = 0.02
// 卡尔曼增益k = P(k | k - 1) / (P(k | k - 1) + R) = 0.0741
// X(k) = 23.9 + 0.0741*(24.1 - 23.9) = 23.915摄氏度。
// k时刻的偏差为P(k) = (1 - K * H)P(k | k - 1) = 0.0186。
// 最后由X(k)和P(k)得出Z(k + 1)。


ArrayXd KF1() {
	int N = 120; // 循环次数
	ArrayXd X(N), Xkf(N), Z(N), P(N);
	X(0) = 25.1; //
	P(0) = 0.01; //
	Z(0) = 24.9; //
	Xkf(0) = Z(0);
	double Q = 0.01; double R = 0.25; // W(k) 的方差 V(k) 的方差
	ArrayXd W = sqrt(Q) * Eigen::ArrayXd::Random(N);
	ArrayXd V = sqrt(R) * Eigen::ArrayXd::Random(N);

	double F = 1.0; double G = 1.0; double H = 1.0; MatrixXd I = MatrixXd::Identity(1, 1);
	for (int kk = 1; kk < N; kk++) {
		X(kk) = F * X(kk - 1) + G * W(kk - 1); // 实际的
		Z(kk) = H * X(kk - 1) + V(kk); // 实际的
		double X_pre = F * Xkf(kk - 1); // 状态预测部分
		double P_pre = F * P(kk - 1) * F + Q; // 协方差矩阵预测
		double Kg = P_pre * 1.0 / (H * P_pre * H + R); // 1x1的元素求逆矩阵 卡尔曼增益更新
		double e = Z(kk) - H * X_pre; // 增益系数更新
		Xkf(kk) = X_pre + Kg * e; // 状态更新
		P(kk) = (I(0, 0) - Kg * H)*P_pre; // 协方差矩阵更新
	}

	return Xkf; // 返回卡尔曼预测
}

// 卡尔曼滤波在自由落体运动目标跟踪的运用
MatrixXd KF2() {
	int N = 1000; // 仿真时间
	MatrixXd Q = MatrixXd::Zero(2, 2); // 过程噪声的方差 = 0
	double R = 1; // 观测噪声方差
	MatrixXd W = sqrt(Q.array()).matrix() * MatrixXd::Random(2, N); // 
	MatrixXd V = sqrt(R) * MatrixXd::Random(2, N); // 观测噪声
	//
	MatrixXd A(2, 2);
	A << 1, 1,
		0, 1; // 状态转移矩阵
	// X ∈ (位移, 速度)
	// X_k+1 = X_k x A
	// 新位移 = 位移 * 1 + 速度 * 1
	// 新速度 = 位移 * 0 + 速度 * 1
	MatrixXd B(2, 1);
	B << 0.5,
		1; // 控制量
	double U = -1; // 控制
	// 位移增加 速度增加
	MatrixXd H(1, 2);
	H << 1, 0; // 观测矩阵
	// -------------------------------------------------------- 初始化 -------------------------------------------------------
	MatrixXd X = MatrixXd::Zero(2, N); // 物体真实运动状态
	X.col(0) << 95, 1; // 初始位移以及速度
	MatrixXd P0(2, 2);
	P0 << 10, 0, 
		0, 1; // 初始误差
	MatrixXd Z = MatrixXd::Zero(1, N);
	Z(0) = (H * X.col(0))(0, 0); // 初始观测值
	MatrixXd Xkf = MatrixXd::Zero(2, N); // 卡尔曼估计状态初始化
	Xkf.col(0) = X.col(0);
	MatrixXd I = MatrixXd::Identity(2, 2);

	for (int kk = 1; kk < N; kk++) {
		X.col(kk) = A * X.col(kk - 1) + B * U; // 真实状态更新
		Z(kk) = (A * P0 * A.transpose() + Q)(0, 0); // 协方差更新
		// 卡尔曼滤波
		MatrixXd X_pre = (A * Xkf.col(kk - 1) + B * U); // 状态预测
		MatrixXd P_pre = (A * P0 * A.transpose() + Q); // 协方差预测
		MatrixXd tmp = (H * P_pre * H.transpose()).array() + R; // 这是一个数字
		MatrixXd Kg = P_pre * H.transpose() * (tmp).inverse(); // 计算卡尔曼增益
		Xkf.col(kk) = X_pre + Kg * (Z(kk) - (H * X_pre)(0,0)); // 状态更新
		P0 = (I - Kg * H) * P_pre; // 方差更新
	}

	return Xkf;
}
