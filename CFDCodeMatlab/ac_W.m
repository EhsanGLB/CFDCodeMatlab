function ac_W = ac_W(I, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fc_w=0.5*rho*u(J,I-1);
Dc_w=D/Ax;
ac_W=Dc_w+max(Fc_w,0);
