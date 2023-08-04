function ac_S = ac_S(I, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fc_s=0.5*rho*v(J-1,I);
Dc_s=D/Ax;
ac_S=Dc_s+max(Fc_s,0);
