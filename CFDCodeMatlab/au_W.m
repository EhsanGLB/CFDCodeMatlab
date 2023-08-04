function au_W = au_W(i, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fu_w=0.5*rho*(u(J,i-1)+u(J,i));
Du_w=mu/Ax;
au_W=Du_w+max(Fu_w,0);
