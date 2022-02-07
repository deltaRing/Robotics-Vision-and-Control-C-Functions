#include "Chapter2.h"
#include <vector>

struct TPoly {
	VectorXd s; // 生成多项式轨迹 
	VectorXd sd; // 速度
	VectorXd sdd; // 加速度

	ArrayXXd S; // 生成多项式轨迹 
	ArrayXXd Sd; // 速度
	ArrayXXd Sdd; // 加速度
};

double ceil(double tmp) {
	if (tmp > 0)
		return tmp - int(tmp) > 0 ? int(tmp) + 1 : int(tmp);
	else if (tmp < 0)
		return int(tmp);
	return 0;
}

ArrayXd ceil(ArrayXd tmp) {
	ArrayXd ttt(tmp.size());
	for (int iii = 0; iii < tmp.size(); iii++) {
		if (tmp(iii) > 0)
			ttt(iii) = tmp(iii) - int(tmp(iii)) > 0 ? int(tmp(iii)) + 1 : int(tmp(iii));
		else if (tmp(iii) < 0)
			ttt(iii) = int(tmp(iii));
		else
			ttt(iii) = 0;
	}
	return ttt;
}

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

// Compute a joint space trajectory
ArrayXXd jtraj(ArrayXd q0, ArrayXd q1, ArrayXd tv, ArrayXd qd0, ArrayXd qd1) {
	ArrayXd t;
	double tscal;
	if (tv.size() > 1) {
		tscal = tv.maxCoeff();
		t = tv / tscal;
	}
	else {
		tscal = 1;
		t = ArrayXd::LinSpaced(int(tv(0)), 0, tv(0) - 1);
	}

	// compute the polynomial coefficients
	ArrayXd A = 6 * (q1 - q0) - 3 * (qd1 + qd0)*tscal;
	ArrayXd B = -15 * (q1 - q0) + (8 * qd0 + 7 * qd1)*tscal;
	ArrayXd C = 10 * (q1 - q0) - (6 * qd0 + 4 * qd1)*tscal;
	ArrayXd E = qd0 * tscal; // as the t vector has been normalized
	ArrayXd F = q0;

	ArrayXXd tt(t.size(), 6);
	tt.col(0) = pow(t, 5);
	tt.col(1) = pow(t, 4);
	tt.col(2) = pow(t, 3);
	tt.col(3) = pow(t, 2);
	tt.col(4) = pow(t, 1);
	tt.col(5) = pow(t, 0);

	ArrayXXd c(6, A.size());
	c.row(0) = A;
	c.row(1) = B;
	c.row(2) = C;
	c.row(3) = ArrayXd::Zero(A.size());
	c.row(4) = E;
	c.row(5) = F;

	ArrayXXd _res_ = tt.matrix() * c.matrix();
	return _res_;
}

// MSTRAJ Multi-segment multi-axis trajectory
ArrayXXd mstraj(ArrayXXd segments, ArrayXd qdmax, ArrayXd tsegment, ArrayXd q0, double dt, ArrayXd Tacc) {
	ArrayXXd R;
	if (q0.size() == 0) {
		q0 = segments.row(0);
		segments = segments.block(1, 0, segments.rows() - 1, segments.cols());
	}
	
	if (qdmax.size() > 0 && tsegment.size() > 0) {
		cout << "Error: Must Specify qdMax or tsegment but not both" << endl;
		return R;
	}

	if (segments.cols() != q0.size()) {
		cout << "Error: WP and Q0 must have same number of columns" << endl;
		return R;
	}

	if (qdmax.size() == 0) {
		if (tsegment.size() != segments.rows()) {
			cout << "Error: Length of TSEG does not match number of segments" << endl;
			return R;
		}
	}

	if (tsegment.size() == 0) {
		if (qdmax.size() == 1) {
			double temp = qdmax(0);
			qdmax = ArrayXd::Ones(segments.cols());
			qdmax = qdmax * temp;
		}
		if (qdmax.size() != segments.cols()) {
			cout << "Length of QDMAX does not match number of axes" << endl;
			return R;
		}
	}

	int ns = segments.rows();
	int nj = segments.cols();

	ArrayXd qd0 = ArrayXd::Zero(nj);
	ArrayXd qdf = ArrayXd::Zero(nj);

	ArrayXd q_prev = q0;
	ArrayXd qd_prev = qd0;
	double clock = 0.0;
	std::vector <double> arrive;
	std::vector <ArrayXd> tg;
	ArrayXd q_next; 
	double tacc2;
	double taccx;
	ArrayXXd qb;

	for (int seg = 0; seg < ns; seg++) {
		double tacc;
		if (Tacc.size() > 1) {
			tacc = Tacc(seg);
		}
		else {
			tacc = Tacc(0);
		}
		tacc = ceil(tacc / dt) * dt;
		tacc2 = ceil(tacc / 2.0 / dt) * dt;
		double taccx;
		if (seg == 0) {
			taccx = tacc2;
		}
		else {
			taccx = tacc;
		}
		
		q_next = segments.row(seg);
		ArrayXd dq = q_next - q_prev;
		ArrayXXd qb;
		double tseg;
		int slowest;

		if (qdmax.size() != 0) {
			qb = taccx * qdmax / 2.0;
			double tb = taccx;

			ArrayXd tl = abs(dq) / qdmax;
			tl = ceil(tl / dt) * dt;

			ArrayXd tt = tb + tl;
			tseg = tt.maxCoeff(&slowest);
			if (tseg <= 2 * tacc)
				tseg = 2 * tacc;
		}
		else if (tsegment.size() > 0) {
			tseg = tsegment(seg);
			slowest = NAN;
		}
		arrive.push_back(clock + tseg);
		if (seg > 0)
			arrive[seg] += tacc2;

		ArrayXd qd = dq / tseg;
		qb = jtraj(q0, q_prev + tacc2 * qd, ArrayXd::LinSpaced(int(taccx / dt), 0, taccx), qd_prev, qd);

		for (int ii = 0; ii < qb.rows() - 1; ii++) {
			tg.push_back(qb.row(ii + 1));
		}
		clock = clock + taccx;
		// add the linear part, from tacc/2+dt to tseg-tacc/2
		for (double t = tacc2 + dt; t < tseg - tacc2; t += dt) {
			double s = t / tseg;
			q0 = (1 - s) * q_prev + s * q_next; // linear step
			tg.push_back(q0);
			clock += dt;
		}
		q_prev = q_next;    // next target becomes previous target
		qd_prev = qd;
	}

	// add the final blend
	qb = jtraj(q0, q_next, ArrayXd::LinSpaced(int(tacc2 / dt), 0, tacc2), qd_prev, qdf);
	for (int ii = 0; ii < qb.rows() - 1; ii++) {
		tg.push_back(qb.row(ii + 1));
	}
	R = ArrayXXd(tg.size(), tg[0].size());

	for (int ii = 0; ii < tg.size(); ii++) {
		R.row(ii) = tg[ii];
	}

	return R;
}

ArrayXd transl(Matrix4d x) {
	// ishomog
	if (abs(x.block(0,0,3,3).determinant() - 1) < 10e-10)
		return x.col(3).segment(0, 3);
	else {
		return Array3d::Zero();
	}
}

ArrayXd tr2delta(Matrix4d A, Matrix4d B) {
	ArrayXd delta;
	Matrix4d TD = A.inverse() * B;

	ArrayXd a1 = transl(TD);
	ArrayXd a2 = vex(t2r(TD) - MatrixXd::Identity(3, 3));

	delta = ArrayXd(a1.size() + a2.size());
	delta.segment(0, a1.size()) = a1;
	delta.segment(a1.size(), a2.size()) = a2;

	return delta;
}
