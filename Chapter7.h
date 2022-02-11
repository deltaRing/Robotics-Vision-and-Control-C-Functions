#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

class Link {
	Link(double theta, double d, double a, double alpha, double sigma);
	Link();
	~Link() {

	}
public:
	double alpha;
	double a;
	double theta;
	double d;
	double sigma;
	double mdh; // standard DH=0, MDH=1
	double offset; // joint coordinate offset
	string name; // joint coordinate name
	bool flip; // joint moves in opposite direction
	MatrixXd qlim; // 2 x 1
	double m; // dynamic: link mass
	MatrixXd r; // 3 x 1
	MatrixXd I; // dynamic: inertia of link with respect to COM (3x3)
	double Jm; // dynamic: motor inertia
	double G; // dynamic: gear ratio
	MatrixXd B; // dynamic: motor viscous friction (1x1 or 2x1)
	MatrixXd Tc; // 2 x 1 dynamic: motor Coulomb friction

public:
	Matrix4d A(double q); // Link transform matrix
};

Link::Link(double theta, double d, double a, double alpha, double sigma) {
	this->theta = theta; this->d = d; this->a = a; this->alpha = alpha; this->sigma = sigma;
	offset = 0;
	flip = false;
	mdh = 0;
	G = 0;
	B = MatrixXd::Zero(1, 1);
	Tc = MatrixXd::Zero(1, 2);
}

Link::Link() {
	alpha = 0;
	a = 0;
	theta = 0;
	d = 0;
	sigma = 0;
	mdh = 0;
	offset = 0;
	flip = false;
	
	m = 0;
	r = MatrixXd::Zero(1, 3);
	I = MatrixXd::Zero(3, 3);
	Jm = 0;
	G = 1;
	B = MatrixXd::Zero(1, 1);
	Tc = MatrixXd::Zero(1, 2);
}

Matrix4d Link::A(double q) {
	double sa = sin(alpha);
	double ca = cos(alpha);
	if (flip)
		q = -q + offset;
	else
		q = q + offset;

	double d, st, ct;
	if (sigma == 0) {
		st = sin(theta);
		ct = cos(theta);
		d = this->d;
	}
	else {
		st = sin(theta);
		ct = cos(theta);
		d = q;
	}

	Matrix4d T;
	if (mdh == 0) {
		T << ct, -st * ca, st * sa, a * ct,
			st, ct * ca, -ct * sa, a * st,
			0, sa, ca, d,
			0, 0, 0, 1;
	}
	else {
		T << ct, -st, 0, a,
			st * ca, ct * ca, -sa, -sa * d,
			st * sa, ct * sa, ca, ca * d,
			0, 0, 0, 1;
	}
	return T;
}

// 连杆坐标系 j-1 -> j 的坐标变换被定义为：
// 基本的旋转和平移
// theta: 关节角度：Z轴绕着X轴的角度旋转
// d: 连杆偏移：z轴的偏移
// a: 连杆长度：沿着xj轴 zj到zj-1轴的距离
// alpha: zj-1轴到z轴之间关于x轴的角度
// sigma: sigma = 0 为转动副，1为移动副
Matrix4d Aj_1_j(double theta, double d, double a, double alpha) {
	Matrix4d _res_;
	_res_ << cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta),
		sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta),
		0, sin(alpha), cos(alpha), d,
		0, 0, 0, 1;
	return _res_;
}
