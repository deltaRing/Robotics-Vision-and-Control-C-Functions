load map1;

free = 1 - map;

free(1,:) = 0;
free(100,:) = 0;
free(:,1) = 0;
free(:,100) = 0;

% skeleton = ithin(free);
% mesh(skeleton) ÎÞ¸Ãº¯Êý
goal = [50,30];        % goal point
start = [20, 10];      % start point
prm = PRM(map);
prm.plan()
prm.plot()
% prm.query(start, goal)