V = diag([0.02, 0.5 * pi / 180] .^ 2); % 这是速度、角度噪声方差
veh = Vehicle(); % 修改vehicle.m 107----109
veh.V = V;
veh.x = [1,1,0];
veh.L = 2; % vehicle.m 没有这个属性需要修改

P0 = diag([0.005, 0.005, 0.001].^2);
ekf = EKF(veh, V, P0);
ekf.run(1000);