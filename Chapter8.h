#include "Chapt2er.h"

MatrixXd tr2jac(Matrix4d T) {
	MatrixXd J = MatrixXd::Zero(6, 6);
	MatrixXd R = t2r(T);
	J.block(0, 0, 3, 3) = R;
	J.block(3, 3, 3, 3) = R;
	return J;
}
// AJB 可以用于计算空间速度 xyz（平移速度） + xyz（旋转速度）


// 计算旋转矩阵导数中的矩阵B，相当于计算空间角速度 P171
Matrix3d rpy2jac(double p, double y, double r, string mode = "xyz") {
	Matrix3d J;
	if (mode == "xyz") {
		J << sin(p), 0, 1,
			-cos(p)*sin(y), cos(y), 0,
			cos(p)*cos(y), sin(y), 0;
	}
	else if (mode == "zyx") {
		J << cos(p)*cos(y), -sin(y), 0,
			cos(p)*cos(y), cos(y),
			0, -sin(p), 0, 1;
	}
	else if (mode == "yxz") {
		J << cos(p)*sin(y), cos(y), 0,
			-sin(p), 0, 1,
			cos(p)*cos(y), -sin(y), 0;
	}
	else {
		J << sin(p), 0, 1,
			-cos(p)*sin(y), cos(y), 0,
			cos(p)*cos(y), sin(y), 0;
	}
	return J;
}

// 上述函数的欧拉角情况
Matrix3d eul2jac(double phi, double theta, double psi) {
	Matrix3d J;
	J << 0, -sin(phi), cos(phi)*sin(theta), 0, cos(phi), sin(phi), sin(theta), 1, 0, cos(theta);
	return J;
}

// MP逆求解
MatrixXd pinv(MatrixXd T) {
	MatrixXd TT;
	JacobiSVD<MatrixXd> svd(T, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();
	MatrixXd A = svd.singularValues();

	int max_size = T.rows() > T.cols() ? T.rows() : T.cols();
	double tol = max_size * 10e-13;
	int r1 = 0;
	for (int ii = 0; ii < A.rows(); ii++) if (A(ii, 0) > tol) r1++;

	MatrixXd _V_(V.rows(), r1);
	MatrixXd _U_(U.rows(), r1);
	MatrixXd _s_(r1, 1);

	for (int ii = 0; ii < r1; ii++) {
		_V_.col(ii) = V.col(ii);
		_U_.col(ii) = U.col(ii);
		_s_(ii, 0) = A(ii, 0);
	}

	MatrixXd __s__ = 1.0 / _s_.array();
	TT = MatrixXd(_V_.rows(), __s__.rows());
	for (int ii = 0; ii < __s__.rows(); ii++) {
		TT.col(ii) = (_V_.array() * __s__(ii, 0)).matrix();
	}
	TT = TT * _U_.transpose();

	return TT;
}
