#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

using namespace Eigen;
using namespace std;

struct AngVec {
	double Theta;
	VectorXd v;
};

struct TrLog {
	double o1;
	VectorXd o2;
};

struct Tr2Rt {
	MatrixXd R;
	VectorXd t;
};

Vector2d RotateAndPosition(double Theta, Vector2d PositionInB, Vector2d Position);
Matrix3d SE2(double x, double y, double Theta); // 相对位姿求解
Matrix3d rotx(double Theta);
Matrix3d roty(double Theta);
Matrix3d rotz(double Theta);
Vector3d tr2eul(Matrix3d R);
Vector3d tr2rpy(Matrix3d R, string order = string("xyz"));
Matrix3d oa2r(Vector3d o, Vector3d a);
Matrix3d rpy2r(double roll, double pitch, double yaw, string order = string("xyz"));
bool isrot(Matrix3d R); // 这个函数是用于判断是否是3x3的矩阵 原作者的代码
bool ishomog(MatrixXd tr); // ISHOMOG Test if SE(3) homogeneous transformation matrix return true if tr = 4x4
MatrixXd t2r(MatrixXd R); // 不知道干啥的 似乎是把4x4 -> 3x3 3x3->2x2
MatrixXd skew(VectorXd v); // skew-symmetric matrix
VectorXd vex(MatrixXd S);
double trace(MatrixXd R);
Tr2Rt tr2rt(MatrixXd T); // Convert homogeneous transform to rotation and translation 
TrLog trlog(MatrixXd T); // TRLOG Logarithm of SO(3) or SE(3) matrix
AngVec tr2angvec(Matrix3d R);
Matrix3d angvec2r(double Theta, Vector3d v);
