function ap_P = ap_P(I, J, u, v, Ax, mu, D)
ap_P=ap_W(I, J, u, v, Ax, mu, D)+ap_E(I, J, u, v, Ax, mu, D)+ap_S(I, J, u, v, Ax, mu, D)+ap_N(I, J, u, v, Ax, mu, D);