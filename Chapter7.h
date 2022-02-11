#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

// 连杆坐标系 j-1 -> j 的坐标变换被定义为：
// 基本的旋转和平移
// theta: 关节角度：Z轴绕着X轴的角度旋转
// d: 连杆偏移：z轴的偏移
// a: 连杆长度：沿着xj轴 zj到zj-1轴的距离
// alpha: zj-1轴到z轴之间关于x轴的角度
Matrix4d Aj_1_j(double theta, double d, double a, double alpha) {
	Matrix4d _res_;
	_res_ << cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta),
		sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta),
		0, sin(alpha), cos(alpha), d,
		0, 0, 0, 1;
	return _res_;
}
