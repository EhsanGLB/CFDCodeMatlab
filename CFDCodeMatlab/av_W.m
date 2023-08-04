function av_W = av_W(I, j, u, v, Ax, mu, D)
%J=j;
i=I;
rho=1000;
%mu=0.0007;
Fv_w=0.5*rho*(u(j+1,i-1)+u(j,i-1));
Dv_w=mu/Ax;
av_W=Dv_w+max(Fv_w,0);