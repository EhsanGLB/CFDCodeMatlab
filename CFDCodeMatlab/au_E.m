function au_E = au_E(i, J, u, v, Ax, mu, D)
rho=1000;
%mu=0.0007;
Fu_e=0.5*rho*(u(J,i)+u(J,i+1));
Du_e=mu/Ax;
au_E=Du_e+max(0,-Fu_e);