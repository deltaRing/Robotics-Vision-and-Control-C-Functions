T1 = SE2(1, 2, 30 / pi * 180);
T2 = SE2(2, 1, 0);
axis([0 5 0 5])
trplot2(T1, 'frame', '1', 'color', 'b')
hold on
trplot2(T2, 'frame', '2', 'color', 'r')
T3 = T2 * T1;
trplot2(T3, 'frame', '3', 'color', 'g')
T4 = T2 * T1;
trplot2(T2, 'frame', '4', 'color', 'c')
P = [3; 2];
P1 = inv(T1);
P1 = P1(1) * [P; 1];