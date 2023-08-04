function au_N = au_N(i, J, u, v, Ax, mu, D)
I=i;
%j=J;
rho=1000;
%mu=0.0007;
Fu_n=0.5*rho*(v(J,I+1)+v(J,I));
Du_n=mu/Ax;
au_N=Du_n+max(0,-Fu_n);