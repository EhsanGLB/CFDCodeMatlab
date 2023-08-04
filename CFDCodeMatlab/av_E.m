function av_E = av_E(I, j, u, v, Ax, mu, D)
%J=j;
i=I;
rho=1000;
%mu=0.0007;
Fv_e=0.5*rho*(u(j+1,i)+u(j,i));
Dv_e=mu/Ax;
av_E=Dv_e+max(0,-Fv_e);