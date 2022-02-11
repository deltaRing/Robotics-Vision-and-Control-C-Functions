#include "Chapter7.h"

class SerialLink {
	SerialLink(Link * L, int num);
	~SerialLink();
public:
	string name;
	string manuf;
	string comment;
	string ikineType;
	MatrixXd links;
	int n;
	double mdh;
	MatrixXd gravity;
	MatrixXd base;
	MatrixXd tool;
	Link * links = NULL;

public:

};

SerialLink::~SerialLink() {

}

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
