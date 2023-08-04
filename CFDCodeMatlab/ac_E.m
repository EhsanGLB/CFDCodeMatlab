function ac_E = ac_E(I, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fc_e=0.5*rho*u(J,I);
Dc_e=D/Ax;
ac_E=Dc_e+max(0,-Fc_e);
