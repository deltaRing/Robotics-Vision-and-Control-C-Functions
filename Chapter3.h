#include "Chapter2.h"

struct TPoly {
	VectorXd s; // 生成多项式轨迹 
	VectorXd sd; // 速度
	VectorXd sdd; // 加速度

	ArrayXXd S; // 生成多项式轨迹 
	ArrayXXd Sd; // 速度
	ArrayXXd Sdd; // 加速度
};

// p(x) = 3x ^ 2 + 2x + 1 at x = 5, 7, and 9:
// p = [3 2 1];
// polyval(p, [5 7 9])
// returns the value of a polynomial P evaluated at X. P
// 第三章只使用到了前两个参数，因此，不做过多编写
VectorXd polyval(VectorXd p, VectorXd x) {
	VectorXd _res_;
	if (p.size() == 0 || x.size() == 0) {
		return _res_;
	}
	int nc = p.size();
	VectorXd y = VectorXd::Zero(x.size());
	if (nc > 0) { // 上面判断了P是否大于0为何还要判断？
		y = VectorXd::Ones(x.size()) * p(0);
	}
	for (int ii = 1; ii < nc; ii++) {
		y = x.array() * y.array() + p.array()(ii);
	}
	return y;
}


TPoly tpoly(double q0, double qf, double t) {
	TPoly _res_;
	ArrayXd tt = ArrayXd::LinSpaced(int(t), 0, t - 1);
	double tf = tt.maxCoeff();
	MatrixXd X = MatrixXd(6, 6);

	X << 0, 0, 0, 0, 0, 1.0,
		pow(tf, 5), pow(tf, 4), pow(tf, 3), pow(tf, 2), pow(tf, 1), 1.0,
		0, 0, 0, 0, 1.0, 0,
		5 * pow(tf, 4), 4 * pow(tf, 3), 3 * pow(tf, 2), 2 * pow(tf, 1), 1.0, 0,
		0, 0, 0, 2.0, 0, 0,
		20 * pow(tf, 3), 12 * pow(tf, 2), 6 * tf, 2, 0, 0;

	double qd0 = 0;
	double qdf = 0;

	ArrayXXd temp1(6, 1);
	ArrayXXd temp2(1, 5);
	ArrayXXd temp3(1, 4);
	temp1 << q0, qf, qd0, qdf, 0, 0;
	temp2 << 5, 4, 3, 2, 1;
	temp3 << 4, 3, 2, 1;
	MatrixXd coeffs = X.inverse() * temp1.matrix();
	coeffs = coeffs.transpose();
	MatrixXd coeffs_d = coeffs.row(0).segment(0, 5).array() * temp2.row(0);
	MatrixXd coeffs_dd = coeffs.row(0).segment(0, 4).array() * temp3.row(0);

	VectorXd p = polyval(coeffs.row(0), tt);
	VectorXd pd = polyval(coeffs_d.row(0), tt);
	VectorXd pdd = polyval(coeffs_dd.row(0), tt);

	_res_.s = p;
	_res_.sd = pd;
	_res_.sdd = pdd;
	return _res_;
}


TPoly tpoly(double q0, double qf, double t, double qd0 = 0, double qdf = 0) {
	TPoly _res_;
	ArrayXd tt = ArrayXd::LinSpaced(int(t), 0, t - 1);
	double tf = tt.maxCoeff();
	MatrixXd X = MatrixXd(6, 6);

	X << 0, 0, 0, 0, 0, 1.0,
		pow(tf, 5), pow(tf, 4), pow(tf, 3), pow(tf, 2), pow(tf, 1), 1.0,
		0, 0, 0, 0, 1.0, 0,
		5 * pow(tf, 4), 4 * pow(tf, 3), 3 * pow(tf, 2), 2 * pow(tf, 1), 1.0, 0,
		0, 0, 0, 2.0, 0, 0,
		20 * pow(tf, 3), 12 * pow(tf, 2), 6 * tf, 2, 0, 0;

	ArrayXXd temp1(6, 1);
	ArrayXXd temp2(1, 5);
	ArrayXXd temp3(1, 4);
	temp1 << q0, qf, qd0, qdf, 0, 0;
	temp2 << 5, 4, 3, 2, 1;
	temp3 << 4, 3, 2, 1;
	MatrixXd coeffs = X.inverse() * temp1.matrix();
	coeffs = coeffs.transpose();
	MatrixXd coeffs_d = coeffs.row(0).segment(0, 5).array() * temp2.row(0);
	MatrixXd coeffs_dd = coeffs.row(0).segment(0, 4).array() * temp3.row(0);

	VectorXd p = polyval(coeffs.row(0), tt);
	VectorXd pd = polyval(coeffs_d.row(0), tt);
	VectorXd pdd = polyval(coeffs_dd.row(0), tt);

	_res_.s = p;
	_res_.sd = pd;
	_res_.sdd = pdd;
	return _res_;
}

// 直线与抛物线混合轨迹
// 如果只输入三个参数，V = (q1-q0)/tf * 1.5;
TPoly lspb(double q0, double q1, double t) {
	TPoly _res_;
	ArrayXd tt = ArrayXd::LinSpaced(int(t), 0, t - 1);
	double tf = tt.maxCoeff();
	double V = (q1 - q0) / tf * 1.5;

	int sign = q1 - q0;
	if (sign > 0) {
		sign = 1;
	}
	else if (sign < 0) {
		sign = -1;
	}
	else {
		sign = 0;
	}
	V = abs(V) * sign;
	if (abs(V) < abs(q1 - q0) / tf) {
		cout << "V too small" << endl;
		return _res_;
	}
	else if (abs(V) > 2 * abs(q1 - q0) / tf) {
		cout << "V too big" << endl;
		return _res_;
	}

	if (q0 == q1) {
		_res_.s = VectorXd::Ones(tt.size()) * q0;
		_res_.sd = VectorXd::Zero(tt.size());
		_res_.sdd = VectorXd::Zero(tt.size());
		return _res_;
	}
	double tb = (q0 - q1 + V * tf) / V;
	double a = V / tb;
	VectorXd p = VectorXd::Zero(tt.size());
	VectorXd pd = VectorXd::Zero(tt.size());
	VectorXd pdd = VectorXd::Zero(tt.size());

	for (int ii = 0; ii < tt.size(); ii++) {
		double ttt = tt(ii);
		if (ttt <= tb) {
			p(ii) = q0 + a / 2 * pow(ttt, 2);
			pd(ii) = a * ttt;
			pdd(ii) = a;
		}
		else if (ttt <= (tf - tb)) {
			// linear motion
			p(ii) = (q1 + q0 - V * tf) / 2 + V * ttt;
			pd(ii) = V;
			pdd(ii) = 0;
		}
		else {
			// final blend
			p(ii) = q1 - a / 2 * pow(tf, 2) + a * tf * ttt - a / 2 * pow(ttt, 2);
			pd(ii) = a * tf - a * ttt;
			pdd(ii) = -a;
		}
	}

	_res_.s = p;
	_res_.sd = pd;
	_res_.sdd = pdd;

	return _res_;
}

// 直线与抛物线混合轨迹
// 如果只输入三个参数，V = (q1-q0)/tf * 1.5;
TPoly lspb(double q0, double q1, double t, double V) {
	TPoly _res_;
	ArrayXd tt = ArrayXd::LinSpaced(int(t), 0, t - 1);
	double tf = tt.maxCoeff();

	int sign = q1 - q0;
	if (sign > 0) {
		sign = 1;
	}
	else if (sign < 0) {
		sign = -1;
	}
	else {
		sign = 0;
	}
	V = abs(V) * sign;
	if (abs(V) < abs(q1 - q0) / tf) {
		cout << "V too small" << endl;
		return _res_;
	}
	else if (abs(V) > 2 * abs(q1 - q0) / tf) {
		cout << "V too big" << endl;
		return _res_;
	}

	if (q0 == q1) {
		_res_.s = VectorXd::Ones(tt.size()) * q0;
		_res_.sd = VectorXd::Zero(tt.size());
		_res_.sdd = VectorXd::Zero(tt.size());
		return _res_;
	}
	double tb = (q0 - q1 + V * tf) / V;
	double a = V / tb;
	VectorXd p = VectorXd::Zero(tt.size());
	VectorXd pd = VectorXd::Zero(tt.size());
	VectorXd pdd = VectorXd::Zero(tt.size());

	for (int ii = 0; ii < tt.size(); ii++) {
		double ttt = tt(ii);
		if (ttt <= tb) {
			p(ii) = q0 + a / 2 * pow(ttt, 2); 
			pd(ii) = a * ttt;
			pdd(ii) = a;
		}
		else if (ttt <= (tf - tb)) {
			// linear motion
			p(ii) = (q1 + q0 - V * tf) / 2 + V * ttt;
			pd(ii) = V;
			pdd(ii) = 0;
		}
		else {
			// final blend
			p(ii) = q1 - a / 2 * pow(tf, 2) + a * tf * ttt - a / 2 * pow(ttt,2);
			pd(ii) = a * tf - a * ttt;
			pdd(ii) = -a;
		}
	}

	_res_.s = p;
	_res_.sd = pd;
	_res_.sdd = pdd;

	return _res_;
}

TPoly mtraj(TPoly(*func)(double, double, double), VectorXd q0, VectorXd qf, double M) {
	TPoly _res_;
	ArrayXXd s(int(M), q0.size());
	ArrayXXd sd(int(M), q0.size());
	ArrayXXd sdd(int(M), q0.size());

	for (int ii = 0; ii < q0.size(); ii++) {
		TPoly _tmp_ = (*func)(q0(ii), qf(ii), M);
		s.col(ii) = _tmp_.s;
		sd.col(ii) = _tmp_.sd;
		sd.col(ii) = _tmp_.sd;
	}

	_res_.S = s;
	_res_.Sd = sd;
	_res_.Sdd = sdd;
	return _res_;
}
