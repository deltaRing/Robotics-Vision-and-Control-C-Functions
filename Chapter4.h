#include "Chapter3.h"

struct Bike {
	double x = 0;
	double y = 0;
	double theta = 0;

	ArrayXd x_path;
	ArrayXd y_path;
	ArrayXd theta_path;
};

void BikeModel(Bike & Bike, double v, double w, double dt = 0.001, double t = 1, double handbrake = 0, double a = 1, double L = 1) {
	if (v > 5) // 速度限制
		v = 5;
	if (v < -5) 
		v = -5;
	if (w > 10) // 角速度限制
		w = 10;
	if (w < -10)
		w = -10;
	if (a > 5) // 加速度限制
		a = 5;
	if (a < -5)
		a = -5;
	
	w = tan(w) / L * v;
	double vx = v * cos(Bike.theta);
	double vy = v * sin(Bike.theta);
	
	double xx = Bike.x;
	double yy = Bike.y;
	double th = Bike.theta;

	Bike.x_path = ArrayXd(int(t / dt));
	Bike.y_path = ArrayXd(int(t / dt));
	Bike.theta_path = ArrayXd(int(t / dt));

	int ii = 0;

	for (double tt = 0; tt < t; tt += dt) {
		Bike.x_path(ii) = xx;
		Bike.y_path(ii) = yy;
		Bike.theta_path(ii) = th;
		ii++;
		xx += vx * dt;
		yy += vy * dt;
		th += w * dt;
		v += a * dt;

		if (handbrake) {
			w = 0;
			v = 0;
		}

		if (th > 60) th = 60;
		if (th < -60) th = -60;
		if (v > 5) // 速度限制
			v = 5;
		if (v < -5)
			v = -5;
		if (w > 10) // 角速度限制
			w = 10;
		if (w < -10)
			w = -10;

		vx = v * cos(th);
		vy = v * sin(th);
	}
	Bike.x = xx;
	Bike.y = yy;
	Bike.theta = th;
}

void LaneChangeMode() {
	double v = 1;
	double w = 0;
	Bike Bike;
	BikeModel(Bike, v, w);
	cout << Bike.x_path;
}
