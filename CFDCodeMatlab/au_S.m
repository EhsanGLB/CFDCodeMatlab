function au_S = au_S(i, J, u, v, Ax, mu, D)
I=i;
%j=J-1;
rho=1000;
%mu=0.0007;
Fu_s=0.5*rho*(v(J-1,I+1)+v(J-1,I));
Du_s=mu/Ax;
au_S=Du_s+max(Fu_s,0);