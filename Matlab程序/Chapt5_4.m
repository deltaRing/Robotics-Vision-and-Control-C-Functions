goal = [0,0,0];
start = [0,2,0];
veh = Bicycle('steermax', 1.2);
rrt = RRT(veh, 'goal', goal, 'range', 5);
rrt.plan(rrt)             % create navigation tree
rrt.query(start, goal)  % animate path from this start location
rrt.plot(rrt)