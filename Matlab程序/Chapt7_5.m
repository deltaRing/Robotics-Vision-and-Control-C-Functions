L = Link([0 0 1 0])
L.offset = pi / 4;

s = 'Tz(L1) Rz(q1) Ry(q2) Ty(L2) Tz(L3) Ry(q3) Tx(L6) Ty(L4) Tz(L5) Rz(q4) Ry(q5) Rz(q6)'
dh = DHFactor(s)