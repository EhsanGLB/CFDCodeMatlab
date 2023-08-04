clc
clear
close all

Ax=0.05e-3;%delta x
Ay=Ax;%delta y
Axc=0.05e-3;%for concentration

rho=1000;%density
mu0=0.0007;%viscosity
mu_inf=0.0007;%
D=3e-8;%diffusion
D0=3e-8;%diffusion
D_inf=3e-8;

F_i1=0.9/(6e10);%flow rate in inlet 1
F_i2=1/(6e10);%flow rate in inlet 2

C_i1=1200;%concentration in inlet 1
C_i2=0;%concentration in inlet 2
Su0=1;%Drug consume rate

w1=0.5e-3;
n_w1=w1/Ax;

h1=0.9e-3;
n_h1=h1/Ax;

h2=1e-3;
n_h2=h2/Ax;

w2=0.55e-3;
n_w2=w2/Ax;

LL=2e-3;
n_l=LL/Ax;

v_i1=F_i1/w1;%Inlet Velocity 1
v_i2=F_i2/w1;%Inlet Velocity 2

u_out=(F_i1+F_i2)/w2;

%/////initial guess
u0=1e-8;
v0=1e-8;
p0=1e-8;
c0=500;

%v on walls and outlet
v_w=v_i1/2;
v_ed=0;
v_eu=0;
v_o=0;

%u on walls and outlet
u_i1=0;
u_i2=0;
u_s=0;
u_n=0;

%//////////////////////////////////////////////////////////////////////////////////
v=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);
%/////////////////////for viscosity(x,y) to v
kv=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);
for I=1:n_w1+n_l+2 %N=5
    
        for j=1:n_h2+n_w2+n_h1+1 %M=5
            
            if (I>=n_w1+2) && (j<n_h2+1)
                mu=mu_inf;
                vv0=0;
            else
                mu=mu0;
                vv0=v0;
            end
            kv(j,I)=mu;
            v(j,I)=vv0;
        end
    
end

for I=2:n_w1+2-1%BC in ilet1 and inlet2
    v(1,I)=v_i2;%inlet2
    v((n_h2+n_w2+n_h1)+1,I)=-v_i1;%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%BC on bottom and up horizental wall
    v(n_h2+1,I)=0;%down
    v((n_h2+n_w2+n_h1)+1,I)=0;%up
end

for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left wall
    v(j,1)=2*v_w-v(j,2);
end

for j=1:n_h2+1%Ghost point for right-bottom vertical wall
    v(j,n_w1+2)=2*v_ed-v(j,n_w1+2-1);
end

for j=(n_h2+n_w2)+1:(n_h2+n_w2+n_h1)+1%Ghost point for right-up vertical wall
    v(j,n_w1+n_l+2)=2*v_eu-v(j,n_w1+n_l+2-1);
end

for j=n_h2+1:(n_h2+n_w2)+1%Ghost point for outlet
    v(j,n_w1+n_l+2)=2*v_o-v(j,n_w1+n_l+2-1);
end

%//////////////////////////////////////////////////////////////////////////////////%%
u=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);
%/////////////////////for viscosity(x,y) to u
ku=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);

for J=1:n_h2+n_w2+n_h1+2 %M=5
        
        for i=1:n_w1+n_l+1 %N=5
            
            if (i>n_w1+1) && (J<=n_h2+2-1)
                mu=mu_inf;
                uu0=0;
            else
                mu=mu0;
                uu0=u0;
            end
            ku(J,i)=mu;
            u(J,i)=uu0;
        end
        
end

for J=n_h2+2:(n_h2+n_w2)+2-1%BC in outlet
   u(J,n_w1+n_l+1)=v(J-1,n_w1+n_l+2-1)-v(J,n_w1+n_l+2-1)+u(J,n_w1+n_l+1-1);
end

for J=2:(n_h2+n_w2+n_h1)+2-1 % BC on left-wall
   u(J,1)=0;%i=1
end

for J=(n_h2+n_w2)+2:(n_h2+n_w2+n_h1)+2-1%BC on right-up-vertical wall
    u(J,n_w1+n_l+1)=0;
end

for J=2:n_h2+2-1%BC on right-down-vertical wall
    u(J,n_w1+1)=0;
end

for i=1:n_w1+1%Ghost point for inlet1 and inlet2
   u(1,i)=2*u_i2-u(1+1,i);%inlet2
   u((n_h2+n_w2+n_h1)+2,i)=2*u_i1-u((n_h2+n_w2+n_h1)+2-1,i);%inlet1
end

for i=n_w1+1:n_w1+n_l+1%Ghost point for up and down horizental wall
   u(n_h2+2-1,i)=2*u_s-u(n_h2+2,i);%down
   u((n_h2+n_w2+n_h1)+2,i)=2*u_n-u((n_h2+n_w2+n_h1)+2-1,i);%up
end

%///////////////////////////////////////////////////////////////////////////////////////
p=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
%/////////////////////for viscosity(x,y) to p
kp=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);

for I=1:n_w1+n_l+2
        for J=1:n_h2+n_w2+n_h1+2 %M=5
            if (I>=n_w1+2) && (J<=n_h2+2-1)
                mu=mu_inf;
                pp0=0;
             else
                mu=mu0;
                pp0=p0;
             end
             kp(J,I)=mu;
             p(J,I)=pp0;
        end
end

for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall
    p(J,1)=p(J,1+1);
end

for I=2:n_w1+2-1%ghost point for inlet1 and inlet2
    p(1,I)=-p(1+1,I);%inlet2
    p((n_h2+n_w2+n_h1)+2,I)=-p((n_h2+n_w2+n_h1)+2-1,I);%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%ghost point for up and down horizental wall
    p(n_h2+2-1,I)=p(n_h2+2,I);%down
    p((n_h2+n_w2+n_h1)+2,I)=p((n_h2+n_w2+n_h1)+2-1,I);%up
end

for J=2:n_h2+2-1%ghost point for right-down-vertical wall
    p(J,n_w1+2)=p(J,n_w1+2-1);
end

for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right-up-vertical wall
    p(J,n_w1+n_l+2)=p(J,n_w1+n_l+2-1);
end

for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
    p(J,n_w1+n_l+2)=0;%-p(J,n_w1+n_l+2-1);
end

%///////////////////////////////////////////////////////////////////////////////////////
c=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
%/////////////////////for viscosity(x,y) to c
kc=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);

for I=1:n_w1+n_l+2
        for J=1:n_h2+n_w2+n_h1+2 %M=5
            if (I>=n_w1+2) && ((J<=n_h2+1) || (J>=n_h2+n_w2+2))
                D=D_inf;
                cc0=0;
             else
                D=D0;
                cc0=c0;
             end
             kc(J,I)=D;
             c(J,I)=cc0;
        end
end

for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall(dc/dx=0)
    c(J,1)=c(J,1+1);
end

for I=2:n_w1+2-1%ghost point for inlet1 and inlet2 wall
    c(1,I)=2*C_i2-c(1+1,I);%inlet2
    c((n_h2+n_w2+n_h1)+2,I)=2*C_i1-c((n_h2+n_w2+n_h1)+2-1,I);%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%ghost point for up and down horizental wall
    c(n_h2+2-1,I)=c(n_h2+2,I);%down
    c((n_h2+n_w2)+2,I)=c((n_h2+n_w2)+2-1,I);%up
end

for J=2:n_h2+2-1%ghost point for right wall down
    c(J,n_w1+2)=c(J,n_w1+n_l+2-1);
end

for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right wall up
    c(J,n_w1+2)=c(J,n_w1+2-1);
end

for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
    c(J,n_w1+n_l+2)=c(J,n_w1+n_l+2-1);
end

%///////////////////////////////////////////
S_u=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);

for I=1:n_w1+n_l+2
        for J=1:n_h2+n_w2+n_h1+2 %M=5
            if (I>=n_w1+2) && ((J<=n_h2+n_w2+1) && (J>=n_h2+2))
                Su=-Su0;
             else
                Su=0;
             end
             S_u(J,I)=Su;
        end
end
%/////////////////////////////////////////////////////////////////////////////////////matrix
uu=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);
erorr_u=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);
r_u=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);

vv=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);
erorr_v=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);
r_v=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);

pp=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
erorr_p=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
b=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);

cc=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
erorr_c=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
r_c=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);

%under relaxation
alpha_u=0.1;
alpha_v=0.1;
alpha_p=0.001;
alpha_c=0.001;

%setting
n=1;
m=5000;
tolorance=1e-3;
Tol=0.1;
mm=2000;
nn=1;

for o=1:m%Tol>tolorance
    
    %/////////////////////////////////////////////////////////////////////////////////////////////////////Step1:Solve discretised momentum equation for u.
    
    for J=2:n_h2+n_w2+n_h1+2-1 %M=5
        
        for i=2:n_w1+n_l+1-1 %N=5
            
            if (i>n_w1+1) && (J<=n_h2+2-1)
                uu(J,i)=0;
            else
                I=i;
                uu(J,i)=(au_W(i,J,u,v,Ax, ku(J,i), D)*u(J,i-1)+au_E(i,J,u,v,Ax, ku(J,i), D)*u(J,i+1)+au_S(i,J,u,v,Ax, ku(J,i), D)*u(J-1,i)+au_N(i,J,u,v,Ax, ku(J,i), D)*u(J+1,i)-(p(J,I+1)-p(J,I)))/(au_P(i,J,u,v,Ax, ku(J,i), D));
            end
            
        end
        
    end
    
    %boundaaey condition of numirical in outlet du/dx=0
    for J=n_h2+2:(n_h2+n_w2)+2-1%BC in outlet
   uu(J,n_w1+n_l+1)=v(J-1,n_w1+n_l+2-1)-v(J,n_w1+n_l+2-1)+uu(J,n_w1+n_l+1-1);
end

for J=2:(n_h2+n_w2+n_h1)+2-1 % BC on left-wall
   uu(J,1)=0;%i=1
end

for J=(n_h2+n_w2)+2:(n_h2+n_w2+n_h1)+2-1%BC on right-up-vertical wall
    uu(J,n_w1+n_l+1)=0;
end

for J=2:n_h2+2-1%BC on right-down-vertical wall
    uu(J,n_w1+1)=0;
end

for i=1:n_w1+1%Ghost point for inlet1 and inlet2
   uu(1,i)=2*u_i2-uu(1+1,i);%inlet2
   uu((n_h2+n_w2+n_h1)+2,i)=2*u_i1-uu((n_h2+n_w2+n_h1)+2-1,i);%inlet1
end

for i=n_w1+1:n_w1+n_l+1%Ghost point for up and down horizental wall
   uu(n_h2+2-1,i)=2*u_s-uu(n_h2+2,i);%down
   uu((n_h2+n_w2+n_h1)+2,i)=2*u_n-uu((n_h2+n_w2+n_h1)+2-1,i);%up
end
    
    %to calculate the erorr_u
    for J=1:n_h2+n_w2+n_h1+2
        for i=1:n_w1+n_l+1
        erorr_u(J,i)=(uu(J,i)-u(J,i))/uu(J,i);
        end
    end
    
    %to calculate the Tol_u
    Tol_u=max(abs(erorr_u));
        
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////Step1:Solve discretised momentum equation for v.
    
    for I=2:n_w1+n_l+2-1 %N=5
    
        for j=2:n_h2+n_w2+n_h1+1-1 %M=5

            if (I>=n_w1+2) && (j<n_h2+1)
                vv(j,I)=0;
            else
                J=j;
                vv(j,I)=(av_W(I,j,u,v,Ax, kv(j,I), D)*v(j,I-1)+av_E(I,j,u,v,Ax, kv(j,I), D)*v(j,I+1)+av_S(I,j,u,v,Ax, kv(j,I), D)*v(j-1,I)+av_N(I,j,u,v,Ax, kv(j,I), D)*v(j+1,I)-(p(J+1,I)-p(J,I)))/(av_P(I,j,u,v,Ax, kv(j,I), D));
            end
            
        end
    
    end
    
    %boundaaey condition in inlets
    for I=2:n_w1+2-1%BC in ilet1 and inlet2
    vv(1,I)=v_i2;%inlet2
    vv((n_h2+n_w2+n_h1)+1,I)=-v_i1;%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%BC on bottom and up horizental wall
    vv(n_h2+1,I)=0;%down
    vv((n_h2+n_w2+n_h1)+1,I)=0;%up
end

for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left wall
    vv(j,1)=2*v_w-vv(j,2);
end

for j=1:n_h2+1%Ghost point for right-bottom vertical wall
    vv(j,n_w1+2)=2*v_ed-vv(j,n_w1+2-1);
end

for j=(n_h2+n_w2)+1:(n_h2+n_w2+n_h1)+1%Ghost point for right-up vertical wall
    vv(j,n_w1+n_l+2)=2*v_eu-vv(j,n_w1+n_l+2-1);
end

for j=n_h2+1:(n_h2+n_w2)+1%Ghost point for outlet
    vv(j,n_w1+n_l+2)=2*v_o-vv(j,n_w1+n_l+2-1);
end
    
    %to calculate the erorr_v
    for I=1:n_w1+n_l+2
        for j=1:n_h2+n_w2+n_h1+1
        erorr_v(j,I)=(vv(j,I)-v(j,I))/vv(j,I);
        end
    end
    
    % to calculate the Tol_v
    Tol_v=max(abs(erorr_v));
    
    %/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Step2: Solve pressure correction equation.
    
    pp(2,2)=(ap_E(2, 2, u, v, Ax, kp(2,2), D)*p(2,2+1)+ap_N(2, 2, u, v, Ax, kp(2,2), D)*p(2+1,2)+b_p(2, 2, u, v, Ax, kp(2,2), D))/(ap_E(2, 2, u, v, Ax, kp(2,2), D)+ap_N(2, 2, u, v, Ax, kp(2,2), D));
        
    pp(n_h2+n_w2+n_h1+2-1,2)=(ap_E(2, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(2,n_h2+n_w2+n_h1+2-1), D)*p(n_h2+n_w2+n_h1+2-1,2+1)+ap_S(2, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(2,n_h2+n_w2+n_h1+2-1), D)*p(n_h2+n_w2+n_h1+2-1,2)+b_p(2, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(2,n_h2+n_w2+n_h1+2-1), D))/(ap_E(2, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(2,n_h2+n_w2+n_h1+2-1), D)+ap_S(2, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(2,n_h2+n_w2+n_h1+2-1), D));
        
    for I=3:n_w1+n_l+2-2 %N=5
    
        pp(2,I)=(ap_W(I, 2, u, v, Ax, kp(2,I), D)*p(2,I-1)+ap_E(I, 2, u, v, Ax, kp(2,I), D)*p(2,I+1)+ap_N(I, 2, u, v, Ax, kp(2,I), D)*p(2+1,I)+b_p(I, 2, u, v, Ax, kp(2,I), D))/(ap_W(I, 2, u, v, Ax, kp(2,I), D)+ap_E(I, 2, u, v, Ax, kp(2,I), D)+ap_N(I, 2, u, v, Ax, kp(2,I), D));
        
        for J=3:n_h2+n_w2+n_h1+2-2 %M=5
        
            pp(J,2)=(ap_E(2, J, u, v, Ax, kp(J,2), D)*p(J,2+1)+ap_S(2, J, u, v, Ax, kp(J,2), D)*p(J-1,2)+ap_N(2, J, u, v, Ax, kp(J,2), D)*p(J+1,2)+b_p(2, J, u, v, Ax, kp(J,2), D))/(ap_E(2, J, u, v, Ax, kp(J,2), D)+ap_S(2, J, u, v, Ax, kp(J,2), D)+ap_N(2, J, u, v, Ax, kp(J,2), D));
            
            pp(J,I)=(ap_W(I, J, u, v, Ax, kp(J,I), D)*p(J,I-1)+ap_E(I, J, u, v, Ax, kp(J,I), D)*p(J,I+1)+ap_S(I, J, u, v, Ax, kp(J,I), D)*p(J-1,I)+ap_N(I, J, u, v, Ax, kp(J,I), D)*p(J+1,I)+b_p(I, J, u, v, Ax, kp(J,I), D))/(ap_P(I, J, u, v, Ax, kp(J,I), D));
    
            pp(J,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D)*p(J,n_w1+n_l+2-1-1)+ap_S(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D)*p(J-1,n_w1+n_l+2-1)+ap_N(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D)*p(J+1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D))/(ap_W(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D)+ap_S(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D)+ap_N(n_w1+n_l+2-1, J, u, v, Ax, kp(J,n_w1+n_l+2-1), D));
    
        end
    
        pp(n_h2+n_w2+n_h1+2-1,I)=(ap_W(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D)*p(n_h2+n_w2+n_h1+2-1,I-1)+ap_E(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D)*p(n_h2+n_w2+n_h1+2-1,I+1)+ap_S(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D)*p(n_h2+n_w2+n_h1+2-1,I)+b_p(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D))/(ap_W(I,n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D)+ap_E(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D)+ap_S(I, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,I), D));
        
    end
    
    pp(2,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, 2, u, v, Ax, kp(2,n_w1+n_l+2-1), D)*p(2,n_w1+n_l+2-1-1)+ap_N(n_w1+n_l+2-1, 2, u, v, Ax, kp(2,n_w1+n_l+2-1), D)*p(2+1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, 2, u, v, Ax, kp(2,n_w1+n_l+2-1), D))/(ap_W(n_w1+n_l+2-1, 2, u, v, Ax, kp(2,n_w1+n_l+2-1), D)+ap_N(n_w1+n_l+2-1, 2, u, v, Ax, kp(2,n_w1+n_l+2-1), D));
    
    pp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1), D)*p(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1-1)+ap_S(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1), D)*p(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1), D))/(ap_W(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1), D)+ap_S(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax, kp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1), D));
     
    
    for I=1:n_w1+n_l+2
        for J=1:n_h2+n_w2+n_h1+2 %M=5
            if (I>=n_w1+2) && (J<=n_h2+1)
                mu=mu_inf;
                pp(J,I)=0;
             else
                mu=mu0;
                pp0=p0;
             end
             kp(J,I)=mu;
        end
    end

    for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall
    pp(J,1)=pp(J,1+1);
end

for I=2:n_w1+2-1%ghost point for inlet1 and inlet2
    pp(1,I)=-pp(1+1,I);%inlet2
    pp((n_h2+n_w2+n_h1)+2,I)=-pp((n_h2+n_w2+n_h1)+2-1,I);%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%ghost point for up and down horizental wall
    pp(n_h2+2-1,I)=pp(n_h2+2,I);%down
    pp((n_h2+n_w2+n_h1)+2,I)=pp((n_h2+n_w2+n_h1)+2-1,I);%up
end

for J=2:n_h2+2-1%ghost point for right-down-vertical wall
    pp(J,n_w1+2)=pp(J,n_w1+2-1);
end

for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right-up-vertical wall
    pp(J,n_w1+n_l+2)=pp(J,n_w1+n_l+2-1);
end

for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
    pp(J,n_w1+n_l+2)=0;%-p(J,n_w1+n_l+2-1);
end

    %to calculate the erorr_p
    for J=1:n_h2+n_w2+n_h1+2
        for I=1:n_w1+n_l+2
        erorr_p(J,I)=(pp(J,I)-p(J,I))/pp(J,I);
        end
    end
    
    % to calculate the Residual_p and Tol_p
    Tol_p=max(abs(erorr_p));
        
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////Step3: Correct pressure and velocities
    
    %this function full the p=p*+p'
    for J=2:n_h2+n_w2+n_h1+2-1
       for I=2:n_w1+n_l+2-1
           p(J,I)=p(J,I)+alpha_p*pp(J,I);%p=p*+alpha_p*p'
       end
    end
    
    for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall
    p(J,1)=p(J,1+1);
end

for I=2:n_w1+2-1%ghost point for inlet1 and inlet2
    p(1,I)=-p(1+1,I);%inlet2
    p((n_h2+n_w2+n_h1)+2,I)=-p((n_h2+n_w2+n_h1)+2-1,I);%inlet1
end

for I=n_w1+2:n_w1+n_l+2-1%ghost point for up and down horizental wall
    p(n_h2+2-1,I)=p(n_h2+2,I);%down
    p((n_h2+n_w2+n_h1)+2,I)=p((n_h2+n_w2+n_h1)+2-1,I);%up
end

for J=2:n_h2+2-1%ghost point for right-down-vertical wall
    p(J,n_w1+2)=p(J,n_w1+2-1);
end

for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right-up-vertical wall
    p(J,n_w1+n_l+2)=p(J,n_w1+n_l+2-1);
end

for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
    p(J,n_w1+n_l+2)=0;%-p(J,n_w1+n_l+2-1);
end

    %this function full the v=alpha_v*(v*+d(p'-p"))+(1-alpha_v)*v matrix on cell surface.
    for I=2:n_w1+n_l+2-1
       for j=2:n_h2+n_w2+n_h1+1-1
           if (I>=n_w1+2) && (j<n_h2+1)
                v(j,I)=0;
            else
                J=j;
                v(j,I)=alpha_v*(vv(j,I)+(pp(J-1,I)-pp(J,I))/av_P(I, j, u, v, Ax, kv(j,I), D))+(1-alpha_v)*v(j,I);
            end
           
       end
    end
    
    for I=2:n_w1+2-1%BC in ilet1 and inlet2
    v(1,I)=v_i2;%inlet2
    v((n_h2+n_w2+n_h1)+1,I)=-v_i1;%inlet1
    end

for I=n_w1+2:n_w1+n_l+2-1%BC on bottom and up horizental wall
    v(n_h2+1,I)=0;%down
    v((n_h2+n_w2+n_h1)+1,I)=0;%up
end

for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left wall
    v(j,1)=2*v_w-v(j,2);
end

for j=1:n_h2+1%Ghost point for right-bottom vertical wall
    v(j,n_w1+2)=2*v_ed-v(j,n_w1+2-1);
end

for j=(n_h2+n_w2)+1:(n_h2+n_w2+n_h1)+1%Ghost point for right-up vertical wall
    v(j,n_w1+n_l+2)=2*v_eu-v(j,n_w1+n_l+2-1);
end

for j=n_h2+1:(n_h2+n_w2)+1%Ghost point for outlet
    v(j,n_w1+n_l+2)=2*v_o-v(j,n_w1+n_l+2-1);
end

    
    %this function full the u=alpha_u*(u*+d(p'-p"))+(1-alpha_u)*u matrix on cell surface.
    for J=2:n_h2+n_w2+n_h1+2-1
       for i=2:n_w1+n_l+1-1
           if (i>n_w1+1) && (J<=n_h2+2-1)
                u(J,i)=0;
            else
                I=i;
                u(J,i)=alpha_u*(uu(J,i)+(pp(J,I-1)-pp(J,I))/au_P(i, J, u, v, Ax, ku(J,i), D))+(1-alpha_u)*u(J,i);
            end
       end
    end
    
    for J=n_h2+2:(n_h2+n_w2)+2-1%BC in outlet
   u(J,n_w1+n_l+1)=v(J-1,n_w1+n_l+2-1)-v(J,n_w1+n_l+2-1)+u(J,n_w1+n_l+1-1);
end

for J=2:(n_h2+n_w2+n_h1)+2-1 % BC on left-wall
   u(J,1)=0;%i=1
end

for J=(n_h2+n_w2)+2:(n_h2+n_w2+n_h1)+2-1%BC on right-up-vertical wall
    u(J,n_w1+n_l+1)=0;
end

for J=2:n_h2+2-1%BC on right-down-vertical wall
    u(J,n_w1+1)=0;
end

for i=1:n_w1+1%Ghost point for inlet1 and inlet2
   u(1,i)=2*u_i2-u(1+1,i);%inlet2
   u((n_h2+n_w2+n_h1)+2,i)=2*u_i1-u((n_h2+n_w2+n_h1)+2-1,i);%inlet1
end

for i=n_w1+1:n_w1+n_l+1%Ghost point for up and down horizental wall
   u(n_h2+2-1,i)=2*u_s-u(n_h2+2,i);%down
   u((n_h2+n_w2+n_h1)+2,i)=2*u_n-u((n_h2+n_w2+n_h1)+2-1,i);%up
end

    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residuals
    
    %////////////////////////%to calculate the Residual_u
     for J=2:n_h2+n_w2+n_h1+2-1 %M=5
        
        for i=2:n_w1+n_l+1-1 %N=5
            
            I=i;
            r_u(J,i)=au_P(i,J,u,v,Ax, ku(J,i), D)*u(J,i)-(au_W(i,J,u,v,Ax, ku(J,i), D)*u(J,i-1)+au_E(i,J,u,v,Ax, ku(J,i), D)*u(J,i+1)+au_S(i,J,u,v,Ax, ku(J,i), D)*u(J-1,i)+au_N(i,J,u,v,Ax, ku(J,i), D)*u(J+1,i))-(p(J,I)-p(J,I+1));
            
        end
        
    end
    
    Residual_u=max(abs(r_u));
    
    %//////////////////////////////////%to calculate the Residual_v
    for I=2:n_w1+n_l+2-1 %N=5
    
        for j=2:n_h2+n_w2+n_h1+1-1 %M=5
            
            
            J=j;
            r_v(j,I)=av_P(I,j,u,v,Ax, kv(j,I), D)*v(j,I)-(av_W(I,j,u,v,Ax, kv(j,I), D)*v(j,I-1)+av_E(I,j,u,v,Ax, kv(j,I), D)*v(j,I+1)+av_S(I,j,u,v,Ax, kv(j,I), D)*v(j-1,I)+av_N(I,j,u,v,Ax, kv(j,I), D)*v(j+1,I))-(p(J,I)-p(J+1,I));
            
        end
    
    end
    
    Residual_v=max(abs(r_v));
    
    %///////////////////////////%to calculate the Residual_con
    for I=2:n_w1+n_l+2-1
        for J=2:n_h2+n_w2+n_h1+2-1
            i=I;
            j=J;
            rho=1000;
            b(I,J)=rho*(u(J,i-1)-u(J,i)+v(j-1,I)-v(j,I));
        end
    end
    
    Residual_con=max(abs(b));
    
    n=n+1;
    Tol=max([Tol_u Tol_v Tol_p]);
end


while mm>nn
    
    for I=2:n_w1+n_l+2-1 %N=5
    
        for J=2:n_h2+n_w2+n_h1+2-1 %M=5
        
            if (I>=n_w1+2) && ((J<=n_h2+1) || (J>=n_h2+n_w2+2))
                D=D_inf;
                cc(J,I)=0;
             else
                D=D0;
                cc(J,I)=(ac_W(I, J, u, v, Axc, mu, kc(J,I))*c(J,I-1)+ac_E(I, J, u, v, Axc, mu, kc(J,I))*c(J,I+1)+ac_S(I, J, u, v, Axc,  mu, kc(J,I))*c(J-1,I)+ac_N(I, J, u, v, Axc,  mu, kc(J+1,I))*c(J+1,I)+S_u(J,I))/(ac_P(I, J, u, v, Axc,  mu, kc(J,I)));
    
             end
             kc(J,I)=D;
             
        end
    
    end
    
    for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall(dc/dx=0)
        c(J,1)=c(J,1+1);
    end

    for I=2:n_w1+2-1%ghost point for inlet1 and inlet2 wall
        c(1,I)=2*C_i2-c(1+1,I);%inlet2
        c((n_h2+n_w2+n_h1)+2,I)=2*C_i1-c((n_h2+n_w2+n_h1)+2-1,I);%inlet1
    end

    for I=n_w1+2:n_w1+n_l+2-1%ghost point for up and down horizental wall
        c(n_h2+2-1,I)=c(n_h2+2,I);%down
        c((n_h2+n_w2)+2,I)=c((n_h2+n_w2)+2-1,I);%up
    end

    for J=2:n_h2+2-1%ghost point for right wall down
        c(J,n_w1+2)=c(J,n_w1+n_l+2-1);
    end

    for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right wall up
        c(J,n_w1+2)=c(J,n_w1+2-1);
    end

    for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
        c(J,n_w1+n_l+2)=c(J,n_w1+n_l+2-1);
    end
    
    %to calculate the erorr_c
    for J=1:n_h2+n_w2+n_h1+2
        for I=1:n_w1+n_l+2
        erorr_c(J,I)=(cc(J,I)-c(J,I))/cc(J,I);
        end
    end
    
    % to calculate the Residual_p and Tol_p
    Tol_c=max(abs(erorr_c));
     %//////////under relaxation this function full the c=c*+alpha_c*c'
    for J=1:n_h2+n_w2+n_h1+2
       for I=1:n_w1+n_l+2
           if (I>=n_w1+2) && ((J<=n_h2+1) || (J>=n_h2+n_w2+2))
                D=D_inf;
                c(J,I)=0;
             else
                D=D0;
                c(J,I)=c(J,I)+alpha_c*cc(J,I);
           end
             kc(J,I)=D;
           
       end
    end
     %////////////////////////%to calculate the Residual_c
     for I=2:n_w1+n_l+2-1 %N=5
    
        for J=2:n_h2+n_w2+n_h1+2-1 %M=5
            
            r_c(J,I)=ac_P(I,J,u,v,Axc, mu, kc(J,I))*c(J,I)-(ac_W(I, J, u, v, Axc, mu, kc(J,I))*c(J,I-1)+ac_E(I, J, u, v, Axc, mu, kc(J,I))*c(J,I+1)+ac_S(I, J, u, v, Axc,  mu, kc(J,I))*c(J-1,I)+ac_N(I, J, u, v, Axc,  mu, kc(J,I))*c(J+1,I)+S_u(J,I));
            
        end
        
    end
    
    Residual_c=max(abs(r_c));
    
    
    nn=nn+1;
end


cont=zeros(n_h2+n_w2+n_h1+2-2,n_w1+n_l+2-2);
for J=2:n_h2+n_w2+n_h1+2-1
    for I=2:n_w1+n_l+2-1
        cont(J,I)=0.5*sqrt( (u(J,I-1)+u(J,I))^2 + (v(J-1,I)+v(J,I))^2 );
    end
end
[X,Y]=meshgrid(1:50,1:50);
surf(X,Y,cont(1:50,1:50))