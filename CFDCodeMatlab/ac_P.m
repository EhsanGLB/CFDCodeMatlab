function ac_P = ac_P(I, J, u, v, Ax, mu, D)
%I=i;
%j=J;
rho=1000;
Fc_w=0.5*rho*u(J,I-1);
Fc_e=0.5*rho*u(J,I);
Fc_s=0.5*rho*v(J-1,I);
Fc_n=0.5*rho*v(J,I);
ac_P=ac_W(I, J, u, v, Ax, mu, D)+ac_E(I, J, u, v, Ax, mu, D)+(Fc_e-Fc_w)+ac_S(I, J, u, v, Ax, mu, D)+ac_N(I, J, u, v, Ax, mu, D)+(Fc_n-Fc_s);