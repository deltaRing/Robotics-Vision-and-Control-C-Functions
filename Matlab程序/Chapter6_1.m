V = diag([0.02, 0.5 * pi / 180] .^ 2); % 这是速度、角度噪声方差
veh = Vehicle(); % 修改vehicle.m 107----109
veh.V = V;
veh.x = [1,1,0];
veh.L = 2; % vehicle.m 没有这个属性需要修改
veh.step(1, 0.3)
veh.add_driver(RandomPath(10))
veh.run()
  
% 注意：我运行的是第一版的程序，上传的是第二版，如果原版本有 请忽略
