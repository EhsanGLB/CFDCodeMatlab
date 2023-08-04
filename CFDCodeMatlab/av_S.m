function av_S = av_S(I, j, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fv_s=0.5*rho*(v(j,I)+v(j-1,I));
Dv_s=mu/Ax;
av_S=Dv_s+max(Fv_s,0);