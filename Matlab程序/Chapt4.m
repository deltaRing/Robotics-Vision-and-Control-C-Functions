% sl_lanechange
sl_drivepoint
xg = [5 5];
x0 = [8 5 pi / 2];
r = sim('sl_drivepoint');
q = find('yout');
plot(q(:,1), q(:,2))