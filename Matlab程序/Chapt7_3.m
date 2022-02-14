mdl_puma560
p560
p560.plot(qr)
qn
T = p560.fkine(qn)
p560.ikine6s(T)
figure
p560.plot(qn)
figure
p560.plot([2.6486   -3.9270    0.0940    2.5326    0.9743    0.3734])
p560.ikine(T)