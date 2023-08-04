function ac_N = ac_N(I, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fc_n=0.5*rho*v(J,I);
Dc_n=D/Ax;
ac_N=Dc_n+max(0,-Fc_n);
