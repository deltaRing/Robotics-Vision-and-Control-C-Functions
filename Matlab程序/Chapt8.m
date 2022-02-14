mdl_puma560
T0 = p560.fkine(qn)
dq = 1e-6
Tp = p560.fkine(qn + [dq 0 0 0 0 0])
dTdq1 = (Tp - T0) / dq
Tp = p560.fkine(qn + [0 dq 0 0 0 0])
dTdq2 = (Tp - T0) / dq

dRdq1 = dTdq1(1:3,1:3);
R = [0         0         1    
         0         1         0  
        -1         0         0 ];
S = dRdq1 * R'
vex(S)

dRdq2 = dTdq2(1:3,1:3);
R = [0         0         1    
         0         1         0  
        -1         0         0 ];
S = dRdq2 *  inv(R)
vex(S)