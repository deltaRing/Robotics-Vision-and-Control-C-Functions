#include "Chapter2.h"

Vector2d RotateAndPosition(double Theta, Vector2d PositionInB, Vector2d Position) {
	MatrixXd _Theta_(2, 2);
	_Theta_ << cos(Theta), -sin(Theta), sin(Theta), cos(Theta); // 定义的旋转矩阵
	return Theta * PositionInB + Position; // 点P自0坐标系到B坐标系的位移
}

Matrix3d SE2(double x, double y, double Theta) {
	Matrix3d _Theta_;
	_Theta_ << cos(Theta), -sin(Theta), x, sin(Theta), cos(Theta), y, 0, 0, 1;
	return _Theta_;
}

Matrix3d rotx(double Theta) {
	Matrix3d _Theta_;
	_Theta_ << 1, 0, 0, 0, cos(Theta), -sin(Theta), 0, sin(Theta), cos(Theta);
	return _Theta_;
}

Matrix3d roty(double Theta) {
	Matrix3d _Theta_;
	_Theta_ << cos(Theta), 0, sin(Theta), 0, 1, 0, -sin(Theta), 0, cos(Theta);
	return _Theta_;
}

Matrix3d rotz(double Theta) {
	Matrix3d _Theta_;
	_Theta_ << cos(Theta), -sin(Theta), 0, sin(Theta), cos(Theta), 0, 0, 0, 1;
	return _Theta_;
}

Vector3d tr2eul(Matrix3d R) {
	Matrix3d RR = R;
	Vector3d eul = Vector3d::Zero();
	if (abs(R(0, 2) < 1e-20 && abs(R(1, 2)) < 1e-20)) {
		eul(0) = 0;
		double sp = 0; double cp = 1;
		eul(1) = atan2(cp*R(0, 2) + sp * R(1, 2), R(2, 2));
		eul(2) = atan2(-sp * R(0, 0) + cp * R(1, 0), -sp * R(0, 1) + cp * R(1, 1));
	}
	else {
		eul(0) = atan2(R(1, 2), R(0, 2));
		double sp = sin(eul(0));
		double cp = cos(eul(0));
		eul(1) = atan2(cp*R(0, 2) + sp * R(1, 2), R(2, 2));
		eul(2) = atan2(-sp * R(0, 0) + cp * R(1, 0), -sp * R(0, 1) + cp * R(1, 1));
	}
	return eul;
}

Vector3d tr2rpy(Matrix3d R, string order) {
	Vector3d rpy = Vector3d::Zero();

	if (order == std::string("xyz") || order == string("XYZ")) {
		if (abs(abs(R(0, 2)) - 1) < 1e-20) {
			rpy(0) = 0;
			if (R(0, 2) > 0)
				rpy(2) = atan2(R(2, 1), R(1, 1));
			else
				rpy(2) = -atan2(R(1, 0), R(2, 1));
			rpy(1) = asin(R(0, 2));
		}
		else
		{
			rpy(0) = -atan2(R(0, 1), R(0, 0));
			rpy(2) = -atan2(R(1, 2), R(2, 2));

			int k = -1; double kk = -1.0;
			k = R(0, 0) > R(0, 1) ? 0 : 1;
			kk = R(0, 0) > R(0, 1) ? R(0, 0) : R(0, 1);
			k = kk > R(1, 2) ? k : 2;
			kk = kk > R(1, 2) ? kk : R(1, 2);
			k = kk > R(2, 2) ? k : 3;
			kk = kk > R(2, 2) ? kk : R(2, 2);
			switch (k)
			{
			case 0:
				rpy(1) = atan(R(0, 2)*cos(rpy(0)) / R(0, 0));
				break;
			case 1:
				rpy(1) = -atan(R(0, 2)*sin(rpy(0)) / R(0, 1));
				break;
			case 2:
				rpy(1) = -atan(R(0, 2)*sin(rpy(0)) / R(1, 2));
				break;
			case 3:
				rpy(1) = atan(R(0, 2)*cos(rpy(0)) / R(2, 2));
				break;
			default:
				break;
			}
		}
	}
	else if (order == std::string("zyx") || order == string("ZYX")) {
		if (abs(abs(R(0, 2)) - 1) < 1e-20) {
			rpy(0) = 0;
			if (R(0, 2) < 0)
				rpy(2) = -atan2(R(0, 1), R(0, 2));
			else
				rpy(2) = atan2(-R(0, 1), -R(0, 2));
			rpy(1) = -asin(R(2, 0));
		}
		else {
			rpy(0) = atan2(R(2, 1), R(2, 2));
			rpy(2) = atan2(R(1, 0), R(0, 0));

			int k = -1; double kk = -1.0;
			k = R(0, 0) > R(1, 0) ? 0 : 1;
			kk = R(0, 0) > R(1, 0) ? R(0, 0) : R(1, 0);
			k = kk > R(2, 1) ? k : 2;
			kk = kk > R(2, 1) ? kk : R(2, 1);
			k = kk > R(2, 2) ? k : 3;
			kk = kk > R(2, 2) ? kk : R(2, 2);
			switch (k)
			{
			case 0:
				rpy(1) = -atan(R(2, 0)*cos(rpy(2)) / R(0, 0));
				break;
			case 1:
				rpy(1) = -atan(R(2, 0)*sin(rpy(2)) / R(1, 0));
				break;
			case 2:
				rpy(1) = -atan(R(2, 0)*sin(rpy(0)) / R(2, 1));
				break;
			case 3:
				rpy(1) = -atan(R(2, 0)*cos(rpy(0)) / R(2, 2));
				break;
			default:
				break;
			}
		}
	}
	else if (order == std::string("yxz") || order == string("YXZ")) {
		if (abs(abs(R(1, 2)) - 1) < 1e-20) {
			rpy(0) = 0;
			if (R(1, 2) < 0)
				rpy(2) = -atan2(R(2, 0), R(0, 0));
			else
				rpy(2) = atan2(-R(2, 0), -R(2, 1));
			rpy(1) = -asin(R(1, 2));
		}
		else {
			rpy(0) = atan2(R(1, 0), R(1, 1));
			rpy(2) = atan2(R(0, 2), R(2, 2));

			int k = -1; double kk = -1.0;
			k = R(1, 0) > R(1, 1) ? 0 : 1;
			kk = R(1, 0) > R(1, 1) ? R(1, 0) : R(1, 1);
			k = kk > R(0, 2) ? k : 2;
			kk = kk > R(0, 2) ? kk : R(0, 2);
			k = kk > R(2, 2) ? k : 3;
			kk = kk > R(2, 2) ? kk : R(2, 2);
			switch (k)
			{
			case 0:
				rpy(1) = -atan(R(1, 2)*sin(rpy(0)) / R(1, 0));
				break;
			case 1:
				rpy(1) = -atan(R(1, 2)*cos(rpy(0)) / R(1, 1));
				break;
			case 2:
				rpy(1) = -atan(R(1, 2)*sin(rpy(0)) / R(1, 2));
				break;
			case 3:
				rpy(1) = -atan(R(1, 2)*cos(rpy(0)) / R(2, 2));
				break;
			default:
				break;
			}
		}
	}
	else {
		cout << "NOT SUPPORT" << endl;
	}
	return rpy;
}

Matrix3d oa2r(Vector3d o, Vector3d a) {
	Vector3d n = o.cross(a);
	o = a.cross(n);
	Matrix3d R;
	R.col(0) << n(0) / abs(n.sum()), n(1) / abs(n.sum()), n(2) / abs(n.sum());
	R.col(1) << o(0) / abs(o.sum()), o(1) / abs(o.sum()), o(2) / abs(o.sum());
	R.col(2) << a(0) / abs(a.sum()), a(1) / abs(a.sum()), a(2) / abs(a.sum());
	return R;
}
