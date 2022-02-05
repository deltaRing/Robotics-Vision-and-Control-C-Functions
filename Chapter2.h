#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

using namespace Eigen;
using namespace std;

Vector2d RotateAndPosition(double Theta, Vector2d PositionInB, Vector2d Position);
Matrix3d SE2(double x, double y, double Theta); // 相对位姿求解
Matrix3d rotx(double Theta);
Matrix3d roty(double Theta);
Matrix3d rotz(double Theta);
Vector3d tr2eul(Matrix3d R);
Matrix3d oa2r(Vector3d o, Vector3d a);
