#include "Chapt2.h"

Matrix4d trotx(double theta) {
	Matrix4d T = Matrix4d::Zero();
	T.block(0, 0, 3, 3) = rotx(theta);
	T(3, 3) = 1;
	return T;
}

Matrix4d troty(double theta) {
	Matrix4d T = Matrix4d::Zero();
	T.block(0, 0, 3, 3) = roty(theta);
	T(3, 3) = 1;
	return T;
}

Matrix4d trotz(double theta) {
	Matrix4d T = Matrix4d::Zero();
	T.block(0, 0, 3, 3) = rotz(theta);
	T(3, 3) = 1;
	return T;
}

Matrix4d transl(double x, double y, double z) {
	Matrix4d T = Matrix4d::Identity();
	T(0, 3) = x;
	T(1, 3) = y;
	T(2, 3) = z;
	return T;
}
