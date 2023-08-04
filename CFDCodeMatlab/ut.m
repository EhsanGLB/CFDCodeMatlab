clc
clear
close all

Ax=0.05e-3;%delta x
Ay=Ax;%delta y

rho=1000;%density
mu=0.0007;%viscosity
D=3e-8;%diffusion
F_i1=0.9/(6e10);%flow rate in inlet 1
F_i2=1/(6e10);%flow rate in inlet 2
C_i1=1200;%concentration in inlet 1
C_i2=0;%concentration in inlet 2
D_con=1;%Drug consume rate

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

u0=0.000005;
v0=0.00001;

%pressure outlet
p0=0.00001;


u=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);

for J=1:n_h2+n_w2+n_h1+2
    for i=1:n_w1+n_l+1
        u(J,i)=u0;
    end
end

for J=n_h2+2:(n_h2+n_w2)+2-1%BC in outlet
   u(J,n_w1+n_l+1)=u_out;
end

for J=2:(n_h2+n_w2+n_h1)+2-1 % BC on left-wall
   u(J,1)=0;%i=1
end

for J=(n_h2+n_w2)+2:(n_h2+n_w2+n_h1)+2-1%BC on right-up-wall
   u(J,n_w1+n_l+1)=0;
end

for J=2:n_h2+2-1%BC on right-down-wall
   u(J,n_w1+n_l+1)=0;
end

for i=1:n_w1+n_l+1%Ghost point for up and down
   u(1,i)=-u(2,i);%down
   u((n_h2+n_w2+n_h1)+2,i)=-u((n_h2+n_w2+n_h1)+2-1,i);%down
end





v=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);

for j=1:n_h2+n_w2+n_h1+1
    for I=1:n_w1+n_l+2
        v(j,I)=v0;
    end
end

for I=2:n_w2+2-1%BC in ilet1 and inlet2
    v(1,I)=v_i1;%inlet2
    v((n_h2+n_w2+n_h1)+1,I)=-v_i2;%inlet1
end

for I=n_w2+2:n_w2+n_l+2-1%BC on bottom and up horizental wall
    v(1,I)=0;%down
    v((n_h2+n_w2+n_h1)+1,I)=0;%up
end

for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left and right wall
    v(j,1)=-v(j,2);
    v(j,n_w1+n_l+2)=-v(j,n_w1+n_l+2-1);
end



p=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
for J=1:n_h2+n_w2+n_h1+2
    for I=1:n_w1+n_l+2
        p(J,I)=p0;
    end
end

for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall
    p(J,1)=p(J,1+1);
end

for I=2:n_w1+n_l+2-1%ghost point for up and down wall
    p(1,I)=p(1+1,I);%down
    p((n_h2+n_w2+n_h1)+2,I)=p((n_h2+n_w2+n_h1)+2-1,I);%up
end

for J=2:n_h2+2-1%ghost point for right wall down
    p(J,n_w1+n_l+2)=p(J,n_w1+n_l+2-1);
end

for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right wall down
    p(J,n_w1+n_l+2)=p(J,n_w1+n_l+2-1);
end

for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
    p(J,n_w1+n_l+2)=p(J,n_w1+n_l+2-1);
end


uu=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+1);
vv=zeros(n_h2+n_w2+n_h1+1,n_w1+n_l+2);
pp=zeros(n_h2+n_w2+n_h1+2,n_w1+n_l+2);
%under relaxation
alpha_u=0.01;
alpha_v=0.01;
alpha_p=0.001;

n=1;
m=3000;
tolorance=1e-3;
Tol=0.1;

while m>n
    
    %/////////////////////////////////////////////////////////////////////////////////////////////////////Step1:Solve discretised momentum equation for u.
    
    for J=2:n_h2+n_w2+n_h1+2-1 %M=5
        
        for i=2:n_w1+n_l+1-1 %N=5
            
            I=i;
            uu(J,i)=(au_W(i,J,u,v,Ax)*u(J,i-1)+au_E(i,J,u,v,Ax)*u(J,i+1)+au_S(i,J,u,v,Ax)*u(J-1,i)+au_N(i,J,u,v,Ax)*u(J+1,i)-(p(J,I+1)-p(J,I)))/(au_P(i,J,u,v,Ax));
            
        end
        
    end
    
    %%boundaaey condition of numirical in outlet du/dx=0
    for J=n_h2+2:(n_h2+n_w2)+2-1%BC in outlet
        uu(J,n_w1+n_l+1)=0.5*(v(J-1,n_w1+n_l+1)+v(J-1,n_w1+n_l+2))-0.5*(v(J,n_w1+n_l+1)+v(J,n_w1+n_l+2))+uu(J,n_w1+n_l+1-1);
    end

    for J=2:(n_h2+n_w2+n_h1)+2-1 % BC on left-wall
        uu(J,1)=0;%i=1
    end

    for J=(n_h2+n_w2)+2:(n_h2+n_w2+n_h1)+2-1%BC on right-up-wall
        uu(J,n_w1+n_l+1)=0;
    end

    for J=2:n_h2+2-1%BC on right-down-wall
        uu(J,n_w1+n_l+1)=0;
    end

    for i=1:n_w1+n_l+1%Ghost point for up and down
        uu(1,i)=-uu(2,i);%down
        uu((n_h2+n_w2+n_h1)+2,i)=-uu((n_h2+n_w2+n_h1)+2-1,i);%down
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
            
            J=j;
            vv(j,I)=(av_W(I,j,u,v,Ax)*v(j,I-1)+av_E(I,j,u,v,Ax)*v(j,I+1)+av_S(I,j,u,v,Ax)*v(j-1,I)+av_N(I,j,u,v,Ax)*v(j+1,I)-(p(J+1,I)-p(J,I)))/(av_P(I,j,u,v,Ax));
            
        end
    
    end
    
    %%boundaaey condition in inlets
    for I=2:n_w2+2-1%BC in ilet1 and inlet2
        vv(1,I)=v_i1;%inlet2
        vv((n_h2+n_w2+n_h1)+1,I)=-v_i2;%inlet1
    end

    for I=n_w2+2:n_w2+n_l+2-1%BC on bottom and up horizental wall
        vv(1,I)=0;%down
        vv((n_h2+n_w2+n_h1)+1,I)=0;%up
    end

    for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left-wall
        vv(j,1)=-vv(j,2);
    end

    for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for right-wall
        vv(j,n_w1+n_l+2)=-vv(j,n_w1+n_l+2-1);
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
    
    pp(2,2)=(ap_E(2, 2, u, v, Ax)*p(2,2+1)+ap_N(2, 2, u, v, Ax)*p(2+1,2)+b_p(2, 2, u, v, Ax))/(ap_E(2, 2, u, v, Ax)+ap_N(2, 2, u, v, Ax));
        
    pp(n_h2+n_w2+n_h1+2-1,2)=(ap_E(2, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,2+1)+ap_S(2, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,2)+b_p(2, n_h2+n_w2+n_h1+2-1, u, v, Ax))/(ap_E(2, n_h2+n_w2+n_h1+2-1, u, v, Ax)+ap_S(2, n_h2+n_w2+n_h1+2-1, u, v, Ax));
        
    for I=3:n_w1+n_l+2-2 %N=5
    
        pp(2,I)=(ap_W(I, 2, u, v, Ax)*p(2,I-1)+ap_E(I, 2, u, v, Ax)*p(2,I+1)+ap_N(I, 2, u, v, Ax)*p(2+1,I)+b_p(I, 2, u, v, Ax))/(ap_W(I, 2, u, v, Ax)+ap_E(I, 2, u, v, Ax)+ap_N(I, 2, u, v, Ax));
        
        for J=3:n_h2+n_w2+n_h1+2-2 %M=5
            
            pp(J,2)=(ap_E(2, J, u, v, Ax)*p(J,2+1)+ap_S(2, J, u, v, Ax)*p(J-1,2)+ap_N(2, J, u, v, Ax)*p(J+1,2)+b_p(2, J, u, v, Ax))/(ap_E(2, J, u, v, Ax)+ap_S(2, J, u, v, Ax)+ap_N(2, J, u, v, Ax));
            
            pp(J,I)=(ap_W(I, J, u, v, Ax)*p(J,I-1)+ap_E(I, J, u, v, Ax)*p(J,I+1)+ap_S(I, J, u, v, Ax)*p(J-1,I)+ap_N(I, J, u, v, Ax)*p(J+1,I)+b_p(I, J, u, v, Ax))/(ap_P(I, J, u, v, Ax));
    
            pp(J,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, J, u, v, Ax)*p(J,n_w1+n_l+2-1-1)+ap_S(n_w1+n_l+2-1, J, u, v, Ax)*p(J-1,n_w1+n_l+2-1)+ap_N(n_w1+n_l+2-1, J, u, v, Ax)*p(J+1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, J, u, v, Ax))/(ap_W(n_w1+n_l+2-1, J, u, v, Ax)+ap_S(n_w1+n_l+2-1, J, u, v, Ax)+ap_N(n_w1+n_l+2-1, J, u, v, Ax));
    
        end
    
        pp(n_h2+n_w2+n_h1+2-1,I)=(ap_W(I, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,I-1)+ap_E(I, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,I+1)+ap_S(I, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,I)+b_p(I, n_h2+n_w2+n_h1+2-1, u, v, Ax))/(ap_W(I,n_h2+n_w2+n_h1+2-1, u, v, Ax)+ap_E(I, n_h2+n_w2+n_h1+2-1, u, v, Ax)+ap_S(I, n_h2+n_w2+n_h1+2-1, u, v, Ax));
        
    end
    
    pp(2,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, 2, u, v, Ax)*p(2,n_w1+n_l+2-1-1)+ap_N(n_w1+n_l+2-1, 2, u, v, Ax)*p(2+1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, 2, u, v, Ax))/(ap_W(n_w1+n_l+2-1, 2, u, v, Ax)+ap_N(n_w1+n_l+2-1, 2, u, v, Ax));
    
    pp(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1)=(ap_W(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1-1)+ap_S(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax)*p(n_h2+n_w2+n_h1+2-1,n_w1+n_l+2-1)+b_p(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax))/(ap_W(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax)+ap_S(n_w1+n_l+2-1, n_h2+n_w2+n_h1+2-1, u, v, Ax));
     
    for J=2:(n_h2+n_w2+n_h1)+2-1%ghost point for left wall
        pp(J,1)=pp(J,1+1);
    end

    for I=2:n_w1+n_l+2-1%ghost point for up and down wall
        pp(1,I)=pp(1+1,I);%down
        pp((n_h2+n_w2+n_h1)+2,I)=pp((n_h2+n_w2+n_h1)+2-1,I);%up
    end

    for J=2:n_h2+2-1%ghost point for right wall down
        pp(J,n_w1+n_l+2)=pp(J,n_w1+n_l+2-1);
    end

    for J=n_h2+n_w2+2:n_h2+n_w2+n_h1+2-1%ghost point for right wall down
        pp(J,n_w1+n_l+2)=pp(J,n_w1+n_l+2-1);
    end

    for J=n_h2+2:n_h2+n_w2+2-1%ghost point for outlet
        pp(J,n_w1+n_l+2)=pp(J,n_w1+n_l+2-1);
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
    for J=1:n_h2+n_w2+n_h1+2
       for I=1:n_w1+n_l+2
           p(J,I)=p(J,I)+alpha_p*pp(J,I);%p=p*+alpha_p*p'
       end
    end
    
    %this function full the u=alpha_u*(u*+d(p'-p"))+(1-alpha_u)*u matrix on cell surface.
    for J=2:n_h2+n_w2+n_h1+2-1
       for i=2:n_w1+n_l+1-1
           I=i;
           u(J,i)=alpha_u*(uu(J,i)+(pp(J,I-1)-pp(J,I))/au_P(i, J, u, v, Ax))+(1-alpha_u)*u(J,i);
       end
    end
    
    for i=1:n_w1+n_l+1%Ghost point for up and down
        u(1,i)=-u(2,i);%down
        u((n_h2+n_w2+n_h1)+2,i)=-u((n_h2+n_w2+n_h1)+2-1,i);%down
    end

    
    
    %this function full the v=alpha_v*(v*+d(p'-p"))+(1-alpha_v)*v matrix on cell surface.
    for I=2:n_w1+n_l+2-1
       for j=2:n_h2+n_w2+n_h1+1-1
           J=j;
           v(j,I)=alpha_v*(vv(j,I)+(pp(J-1,I)-pp(J,I))/av_P(I, j, u, v, Ax))+(1-alpha_v)*v(j,I);
       end
    end
    
    for j=1:(n_h2+n_w2+n_h1)+1%Ghost point for left and right wall
        v(j,1)=-v(j,2);
        v(j,n_w1+n_l+2)=-v(j,n_w1+n_l+2-1);
    end
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Residuals
    
    %////////////////////////%to calculate the Residual_u
     for J=2:n_h2+n_w2+n_h1+2-1 %M=5
        
        for i=2:n_w1+n_l+1-1 %N=5
            
            I=i;
            r_u(J,i)=au_P(i,J,u,v,Ax)*u(J,i)-(au_W(i,J,u,v,Ax)*u(J,i-1)+au_E(i,J,u,v,Ax)*u(J,i+1)+au_S(i,J,u,v,Ax)*u(J-1,i)+au_N(i,J,u,v,Ax)*u(J+1,i))-(p(J,I)-p(J,I+1));
            
        end
        
    end
    
    Residual_u=max(abs(r_u));
    
    %//////////////////////////////////%to calculate the Residual_v
    for I=2:n_w1+n_l+2-1 %N=5
    
        for j=2:n_h2+n_w2+n_h1+1-1 %M=5
            
            J=j;
            r_v(j,I)=av_P(I,j,u,v,Ax)*v(j,I)-(av_W(I,j,u,v,Ax)*v(j,I-1)+av_E(I,j,u,v,Ax)*v(j,I+1)+av_S(I,j,u,v,Ax)*v(j-1,I)+av_N(I,j,u,v,Ax)*v(j+1,I))-(p(J,I)-p(J+1,I));
            
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
