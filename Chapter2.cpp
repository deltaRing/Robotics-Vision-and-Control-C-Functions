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

Matrix3d rpy2r(double roll, double pitch, double yaw, string order) {
	Matrix3d R;
	if (order == string("xyz") || order == string("XYZ")) {
		R = rotx(yaw) * roty(pitch) * rotz(roll);
	}
	else if (order == string("zyx") || order == string("ZYX")) {
		R = rotz(yaw) * roty(pitch) * rotx(roll);
	}
	else if (order == string("yxz") || order == string("YXZ")) {
		R = roty(yaw) * rotx(pitch) * rotz(roll);
	}
	else {
		cout << "ERROR ORDER" << endl;
	}
	return R;
}

bool isrot(Matrix3d R) {
	if (R.cols() == 3 && R.rows() == 3) {
		/*if (abs(R.determinant() - 1) < 1e-10)
			return true;
		else
			return false;*/ // 可能不对哈
		return true;
	}
	return false;
}

bool ishomog(MatrixXd tr) {
	if (tr.cols() == 4 && tr.rows() == 4)
		return true;
	return false;
}

MatrixXd t2r(MatrixXd T) {
	int d1 = T.rows();
	int d2 = T.cols();
	MatrixXd R;
	if (d1 != d2) {
		cout << "Error: Matrix must be square" << endl;
		return R;
	}
	else if (d1 != 3 || d1 != 4) {
		cout << "Error: Argument is not a homogeneous transform (sequence)" << endl;
		return R;
	}
	int n = d1;
	R = T.block(0, 0, n - 1, n - 1);
	return R;
}

MatrixXd skew(VectorXd v) {
	MatrixXd R;
	if (v.size() == 3) {
		R = MatrixXd(3, 3);
		R << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	}
	else if (v.size() == 1) {
		R = MatrixXd(2, 2);
		R << 0, -v(0), v(0), 0;
	}
	else {
		cout << "Error: SKEW BAD Arguments：V must be 1-v 3-v Vector";
	}
	return R;
}

double trace(MatrixXd R) {
	int ii;
	if (R.cols() != R.rows()) {
		ii = min(R.cols(), R.rows());
	}
	else {
		ii = R.rows();
	}
	double _res_ = 0.0;
	for (int i = 0; i < ii; i++)
		_res_ += R(i, i);
	return _res_;
}

VectorXd vex(MatrixXd S) {
	VectorXd v;
	if (S.rows() == 3 && S.cols() == 3) {
		v = VectorXd(3);
		v << 0.5 * S(2, 1) - 0.5 * S(1, 2), 0.5 * S(0, 2) - 0.5 * S(2, 0), 0.5 * S(1, 0) - 0.5 * S(0, 1);
	}
	else if (S.rows() == 2 && S.cols() == 2) {
		v = VectorXd(1);
		v << 0.5 * S(1, 0) - 0.5 * S(0, 1);
	}
	return v;
}

Tr2Rt tr2rt(MatrixXd T) {
	Tr2Rt _res_;
	if (T.rows() != T.cols()) {
		cout << "Error: T must be a square." << endl;
		return _res_;
	}

	int n = T.rows() - 1;
	_res_.R = T.block(0, 0, n, n);
	_res_.t = T.col(n).segment(0, n);
	return _res_;
}

TrLog trlog(MatrixXd T) {
	TrLog _result_;
	MatrixXd R;
	VectorXd w;
	double theta;
	if (isrot(T)) {
		R = T;
		if (abs(trace(R) - 3) < 10e-20 * 100) {
			w = Vector3d::Zero();
			theta = 0.0;
		}
		else if (abs(trace(R) + 1) < 10e-20 * 100) {
			Vector3d DiagR;
			DiagR << R(0, 0), R(1, 1), R(2, 2);
			int k;
			double mx = DiagR.maxCoeff(&k);
			Matrix3d I = Matrix3d::Identity();

			Vector3d col = R.col(k) + I.col(k);
			w = col / sqrt(2.0 * (1.0 + mx));
			theta = 3.1415926535;
		}
		else {
			theta = acos((trace(R) - 1) / 2.0);
			MatrixXd skw = (R - R.transpose()) / 2.0 / sin(theta);
			w = vex(skw); // is a unit vector
		}
		_result_.o1 = theta;
		_result_.o2 = w;
	}
	else if (ishomog(T)) {
		Tr2Rt Rt = tr2rt(T);
		VectorXd v;
		if (abs(trace(R) - 3) < 10e-20 * 100) {
			w = Vector3d::Zero();
			v = Rt.t;
			theta = 1.0;
			Matrix3d skw = Matrix3d::Zero();
		}
		else {
			TrLog TW = trlog(R);
			theta = TW.o1;
			w = TW.o2;
			MatrixXd skw = skew(w);

			MatrixXd tmp1 = MatrixXd::Identity(skw.rows(), skw.cols()) / theta;
			MatrixXd tmp2 = skw / 2.0;
			MatrixXd tmp3 = (1.0 / theta - tan(3.1415926535 / 2 - (theta / 2)) / 2.0) * skw * skw;

			MatrixXd Ginv = tmp1 - tmp2 + tmp3;
			v = (Ginv * Rt.t).array();
		}
		_result_.o1 = theta;
		_result_.o2 = VectorXd(v.size() + w.size());
		_result_.o2.segment(0, v.size()) = v;
		_result_.o2.segment(v.size(), w.size()) = w;
	}
	else {
		cout << "Error at trlog: Expect SE3 or SO2" << endl;
	}
	return _result_;
}

AngVec tr2angvec(Matrix3d R) {
	AngVec _res_;
	if (!isrot(R)) {
		R = t2r(R);
	}
	if (abs(R.determinant() - 1) > 10e-19) {
		cout << "Error: R is not orthonormal" << endl;
		return _res_;
	}

	TrLog ThV = trlog(R);
	double theta = ThV.o1;
	VectorXd n = ThV.o2;
	_res_.Theta = theta;
	_res_.v = n;
	return _res_;
}

Matrix3d angvec2r(double Theta, Vector3d v) {
	Matrix3d R;
	if (v.norm() < 10e-20) {
		if (abs(Theta) > 0) {
			cout << "Error: norm of direction is zero";
			return R;
		}
		else {
			return Matrix3d::Identity();
		}
	}
	MatrixXd sk = skew(v / v.norm());

	MatrixXd tmp1 = MatrixXd::Identity(R.rows(), R.cols());
	MatrixXd tmp2 = sin(Theta) * sk;
	MatrixXd tmp3 = (1 - cos(Theta)) * sk * sk;

	R = tmp1 + tmp2 + tmp3;
	return R;
}
