function b_p = b_p(I, J, u, v, Ax, mu, D)
i=I;
j=J;
rho=1000;
b_p=rho*(u(J,i-1)-u(J,i)+v(j-1,I)-v(j,I));