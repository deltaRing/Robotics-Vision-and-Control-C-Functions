L(1) = Link([0 0 1 1])
L(2) = Link([0 0 1 1])
two_link = SerialLink(L, 'name', 'two_link')
mdl_planar2
two_link.fkine([pi/4 -pi/4])
two_link.plot([-pi/4 pi/4])