R = rotx(pi / 3) * roty(pi / 3) * rotz(pi / 3);
% trplot(R)
% tranimate(R);
tr2rpy(R);

q = Quaternion(rpy2tr(0.1, 0.2, 0.3))
R = rpy2r(0.1, 0.2, 0.3);
angvec2r(pi/2, [1,0,0]);
tr2angvec(R);

trotx(pi / 2);

a = [1 0 0]';
o = [0 1 0]';
oa2r(o, a);