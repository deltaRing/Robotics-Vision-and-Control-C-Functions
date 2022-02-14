% 若无法运行请再运行一次rvctools中的startup_rtb.m

V = diag([0.02, 0.5 * pi / 180] .^ 2); % 这是速度、角度噪声方差
veh = Vehicle(); % 修改vehicle.m 107----109
veh.V = V;
veh.x = [1;1;0];
veh.L = 2; % vehicle.m 没有这个属性需要修改

map = Map(20, 10);
map.plot()
W = diag([0.1, 1 * pi / 180].^2);
sensor = RangeBearingSensor(veh, map, W);

[z, i] = sensor.reading()
map.feature(12)  % 第十二个障碍物