function av_N = av_N(I, j, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fv_n=0.5*rho*(v(j,I)+v(j+1,I));
Dv_n=mu/Ax;
av_N=Dv_n+max(0,-Fv_n);