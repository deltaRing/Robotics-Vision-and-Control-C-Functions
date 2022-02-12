#include "SE3.h"
#include "Chapter3.h"
#include "Chapter7.h"

class SerialLink {
	SerialLink(Link * L, int num);
public:
	string name;
	string manuf;
	string comment;
	string ikineType;
	int n;
	double mdh;
	MatrixXd gravity;
	MatrixXd base;
	MatrixXd tool;
	MatrixXd qlim;
	Link * links = NULL;

public:
	// bool isspherical(); Link 类 无法使用
	ArrayXi isrevolute() {
		ArrayXi _res_ = ArrayXi::Zero(n);
		if (this->links == NULL)
			return _res_;
		for (int i = 0; i < n; i++) {
			_res_(i) = this->links[i].isrevolute;
		}
		return _res_;
	}
};

SerialLink::SerialLink(Link * L, int num) {
	name = string("noname");
	n = 0;
	mdh = 0;
	gravity = MatrixXd(3, 1);
	gravity << 0, 0, 9.81;
	base = MatrixXd::Identity(4, 4);
	tool = MatrixXd::Identity(4, 4);
	links = L;
	n = num;
	mdh = links[0].mdh;
}

// 这个似乎是给Revolute用的 Link类用不了
/*bool SerialLink::isspherical() {
	Link * L = &this->links[2];

	double pi = 3.1415926535;
	MatrixXd alpha(1, 2);
	alpha << -pi / 2, pi / 2;

	bool tmp1 = L[0].a == 0 && L[1].a == 0;
	bool tmp2 = L[1].d == 0;
	bool tmp3 = L[0].alpha == alpha
}*/

// SerialLink.fkine Forward kinematics
// 正运动学
MatrixXd fkine(SerialLink robot, ArrayXd q) {
	int n = robot.n;
	Link * L = robot.links;
	MatrixXd t = robot.base;
	if (q.size() == n) {
		for (int i = 0; i < n; i++) {
			t = t * L[i].A(q(i));
		}
		t = t * robot.tool;
	}
	return t;
}

// 逆运动学
/*MatrixXd ikine6s(SerialLink robot, Matrix4d TT) {
	MatrixXd thetavec;
	if (robot.mdh != 0) {
		cout << "Error: Solution only applicable for standard DH conventions" << endl;
		return thetavec;
	}

	if (robot.n != 6) {
		cout << "Error: Only support for 6-axis robot" << endl;
		return thetavec;
	}

	TT = SE3(TT);
	Link * L = robot.links;

}*/ // 这是给puma560 Revolute 用的

// 计算雅可比矩阵
MatrixXd jacobe(SerialLink robot, MatrixXd q) {
	MatrixXd J;
	int n = robot.n;
	Link * L = robot.links;

	J = MatrixXd::Zero(6, n);
	MatrixXd U = robot.tool;

	for (int j = n - 1; j >= 0; j--) {
		if (robot.mdh == 0)
			U = L[j].A(q(0, j)) * U;
		MatrixXd UT = U.transpose();
		MatrixXd d(1, 3);
		MatrixXd delta(3, 1);
		if (L[j].isrevolute) {
			d << -UT(0, 0) * UT(1, 3) + UT(1, 0) * UT(0, 3),
				-UT(0, 1) * UT(1, 3) + UT(1, 1) * UT(0, 3),
				-UT(0, 2) * UT(1, 3) + UT(1, 2) * UT(0, 3);
			delta = UT.row(2).segment(0, 3);
		}
		else {
			d = UT.row(2).segment(0, 3);
			delta << 0, 0, 0;
		}
		J.col(j).segment(0, 3) = d;
		J.col(j).segment(3, 3) = delta;

		if (robot.mdh != 0) {
			U = L[j].A(q(0, j)) * U;
		}
	}
	return J;
}

struct ikine_opt {
	double ilimit = 500;
	double rlimit = 100;
	double slimit = 100;
	double tol = 1e-10;
	double lambda = 0.1;
	double lambdamin = 0;
	bool search = false;
	bool quiet = false;
	bool verbose = false;
	MatrixXd mask = MatrixXd::Ones(1, 6);
	MatrixXd q0;
	double transpose = NAN;
};

// 逆运动学
MatrixXd ikine(SerialLink robot, Matrix4d tr) {
	int n = robot.n;
	Matrix4d TT = tr;

	double pi = 3.1415926535;

	ikine_opt opt;
	opt.q0 = MatrixXd::Zero(1, n);
	MatrixXd q(1, n);

	if (opt.search) {
		opt.search = false;
		opt.quiet = true;
		for (int k = 0; k < opt.slimit; k++) {
			for (int j = 0; j < n; j++) {
				MatrixXd qlim = robot.links[j].qlim;
				if (qlim.cols() == 0 || qlim.rows() == 0) {
					if (robot.links[j].isrevolute) {
						q(0, j) = ArrayXd::Random(1)(0) * 2 * pi - pi;
					}
					else {
						cout << "Error: For a prismatic joint, search requires joint limits" << endl;
						return q;
					}
				}
				else {
					q(0, j) = ArrayXd::Random(1)(0) * (qlim(1) - qlim(0)) + qlim(0);
				}
			}


			//
			// q = robot.ikine(tr, q, args{:}, 'setopt', opt); 
			//
		}
		return q; // TODO 
	}

	MatrixXd W = MatrixXd::Identity(opt.mask.size(), opt.mask.size());
	MatrixXd qt = MatrixXd::Zero(1, n);
	int tcount = 0; // total iteration count
	int rejcount = 0; // rejected step count
	q = opt.q0;
	bool failed = false;
	ArrayXi revolutes = robot.isrevolute();

	for (int i = 0; i < 1; i++) {
		Matrix4d T = TT; 
		double lambda = opt.lambda;
		int iterations = 0;

		while (true) {
			iterations += 1;
			if (iterations > opt.ilimit) {
				failed = true;
				break;
			}

			ArrayXd e = tr2delta(fkine(robot, q.array()), T);

			if ((W * e.matrix()).norm() < opt.tol) {
				break;
			}

			// 计算雅可比矩阵
			MatrixXd J = jacobe(robot, q);
			MatrixXd JtJ = J.transpose() * W * J;

			if (!isnan(opt.transpose)) {
				MatrixXd dq = opt.transpose * J.transpose() * e.matrix();
				q = q + dq.transpose();
			}
			else {
				MatrixXd dq = (JtJ + (lambda + opt.lambdamin) * MatrixXd::Identity(JtJ.rows(), JtJ.cols())).inverse() * J.transpose() * W * e.matrix();
				MatrixXd qnew = q + dq.transpose();
				ArrayXd enew = tr2delta(fkine(robot, qnew.array()), T);
			
				if ((W * enew.matrix()).norm() < (W * e.matrix()).norm()) {
					q = qnew;
					e = enew;
					lambda = lambda / 2;
					rejcount = 0;
				}
				else {
					lambda = lambda * 2;
					rejcount = rejcount + 1;
					if (rejcount > opt.rlimit) {
						failed = true;
						break;
					}
					continue;
				}
			}
			for (int kkk = 0; kkk < revolutes.size(); kkk++) {
				if (q(kkk) > pi && revolutes(kkk)) {
					q(kkk) = q(kkk) - 2 * pi;
				}
				if (q(kkk) < -pi && revolutes(kkk)) {
					q(kkk) = q(kkk) + 2 * pi;
				}
			}
		}
	}

	if (failed) {
		cout << "Error: Iteration is failed during IKine" << endl;
		return q;
	}
	return q;
}
