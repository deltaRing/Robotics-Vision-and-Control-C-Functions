% s = tpoly(0, 1, 50)
% plot(s)
% lspb(0, 1, 50, 0.025);
% x = mtraj(@tpoly, [0, 2], [1 -1], 50)
% via = [4,1;4,4;5,2;2,5];
% q = mstraj(via, [2,1], [], [4,1], 0.05, 1);
T0 = transl(1,2,3) * trotx(1) * troty(1) * trotz(1);
T1 = T0 * transl(0.01, 0.02, 0.03) * trotx(0.001) * troty(0.002) * trotz(0.003);
d = tr2delta(T0, T1);
delta2tr(d) * T0