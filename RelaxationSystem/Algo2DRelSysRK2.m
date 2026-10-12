%% Relaxationsys in 2D - Double Shear Layer test case
% Feb/26
% Relaxationssystem mit p Formulierung auf [0,1]^2
% IMEX RK2

clear; close all; clc;
% function Algo2DRelSysRK2(n)
eps = 1e-04; % Geht bis eps = 1e-01

tau = 0.0001;
delta = 0.05;
% rho = pi/15;
% Skaliere die Schichtdicke: pi/15 im [0,2pi]-Gebiet entspricht 1/30 im [0,1]-Gebiet
rho = 1/30;
tic;
tFinal = 0.1%1.5;
CFL = 0.75;

d = 2; % dimension

n = 129; % number of grid pnts


% IMEX Verfahren - Typ A - erste Zeile d. Tableaus 0
c = 1-1/sqrt(2); gamma = (c-1/2)/(c-1);

a21_ = c; 
a31_ = 1-1/(2*c); a32_ = 1/(2*c);

a11 = 0;
a21 = 0; a22 = c;
a31 = 0; a32 = 1-gamma; a33 = gamma;

% Discretize spatial domain
x_min = 0; x_max = 1;
x = linspace(x_min,x_max,n); y = x; dx = x(2)-x(1); dy = dx;

% meshgrid
[X,Y] = meshgrid(x,y); 

%% initial conditions
% initial condition like in Lattice-Boltzmann paper
% u1 = (Y <= pi).*tanh(1/rho*(Y-pi/2)) + (Y> pi).*tanh(1/rho*(3*pi/2-Y));
% Skalierte Initial Condition: 
% pi -> 0.5 (Mitte)
% pi/2 -> 0.25 (untere Scherschicht)
% 3*pi/2 -> 0.75 (obere Scherschicht)
u1 = (Y <= 0.5).*tanh(1/rho*(Y-0.25)) + (Y > 0.5).*tanh(1/rho*(0.75-Y));
% u2 = delta*sin(X);
% Störung u2: Sinuswelle muss über das Gebiet [0,1] eine volle Periode haben (Daher * 2*pi)
u2 = delta*sin(2*pi*X);

% initialize Theta
Theta11 = u1*0; Theta12 = u1*0; Theta21 = u1*0; Theta22 = u1*0;
% initialize Gamma
Gamma1 = u1*0; Gamma2 = u1*0;
% initialize rho, 
rho = ones(n,n); 
% initialize T
T = zeros(n,n);
% initialize p
p = rho + T;

t = 0; dt = dx^2*CFL; 

% Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); A_(1,N-1) = 1; A_(N,2) = 1;
    Lp = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    alpha_e = 1; alpha_i = 0;

    %% --------------- STAGE 1 -------------------------------------------
    u1_1 = u1; u2_1 = u2; rho_1 = rho; p_1 = p; T_1 = T;
    Theta11_1 = Theta11; Theta12_1 = Theta12; Theta21_1 =  Theta21; Theta22_1 = Theta22;
    Gamma1_1 = Gamma1; Gamma2_1 = Gamma2;

    %% --------------- STAGE 2 -------------------------------------------
    %% Compute Gamma^n+1
    [T_dx1,T_dy1] = Grad_O2(T_1,dx,dy);

    Gamma1_2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma1 + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u1_1.*T_1 - T_dx1*tau);
    Gamma2_2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma2 + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u2_1.*T_1 - T_dy1*tau);

    %% Compute p^n+1/2 and T^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_Gamma1_2 = Div_O2(Gamma1_2,Gamma2_2,dx,dy);
    flux_Gamma2_2 = flux_Gamma1_2;
    p_2 = p - dt/d*a22*flux_Gamma2_2/eps; % computation of p^{n+1/2}/eps
    p_2star = p_2; % introduce helper var.
    T_2 = T - dt/d*a22*flux_Gamma1_2; % computation of T^{n+1/2}


    %% Compute Theta^n+1
    Du11_1 = Div_O2(2*u1_1,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_1,0*u1,Theta11,alpha_i,dx,dy);
    Du12_1 = Div_O2(u2_1,u1_1,dx,dy);%Lax_Friedrichs1(u2_1,u1_1,Theta12,alpha_i,dx,dy);
    Du21_1 = Div_O2(u2_1,u1_1,dx,dy);%Lax_Friedrichs1(u2_1,u1_1,Theta21,alpha_i,dx,dy);
    Du22_1 = Div_O2(u2*0,2*u2_1,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_1,Theta22,alpha_i,dx,dy);


    Au11_1 = u1_1.^2 - (u1_1.^2+u2_1.^2)/2;
    Au12_1 = u1_1.*u2_1;
    Au21_1 = u1_1.*u2_1;
    Au22_1 = u2_1.^2 - (u1_1.^2+u2_1.^2)/2;


    div_u_1 = Div_O2(u1_1,u2_1,dx,dy);


    Theta11_2 = eps^2*tau/(eps^2*tau+dt*a22)*Theta11 +...
        +dt/(eps^2*tau+dt*a22)*a21_*(Au11_1 - Du11_1*tau) ...
        + dt*tau/(eps^2*tau+dt*a22)*2/d*a21_*div_u_1;

    Theta12_2 = eps^2*tau/(eps^2*tau+dt*a22)*Theta12 +...
        +dt/(eps^2*tau+dt*a22)*a21_*(Au12_1 - Du12_1*tau);

    Theta21_2 = eps^2*tau/(eps^2*tau+dt*a22)*Theta21 +...
        +dt/(eps^2*tau+dt*a22)*a21_*(Au21_1 - Du21_1*tau);

    Theta22_2 = eps^2*tau/(eps^2*tau+dt*a22)*Theta22 +...
        +dt/(eps^2*tau+dt*a22)*a21_*(Au22_1 - Du22_1*tau) ...
        + dt*tau/(eps^2*tau+dt*a22)*2/d*a21_*div_u_1;


    %% p-Form
    %% Compute u^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_u1_2 = Div_O2(Theta11_2,Theta12_2,dx,dy);%Lax_Friedrichs1(Theta11_2,Theta12_2,u1,alpha_i,dx,dy);
    flux_u2_2 = Div_O2(Theta21_2,Theta22_2,dx,dy);%Lax_Friedrichs1(Theta21_2,Theta22_2,u2,alpha_i,dx,dy);

    u1_2 = u1 - dt*a22*flux_u1_2;
    u2_2 = u2 - dt*a22*flux_u2_2;


    %% Compute p^n+1   
    Lp2 = (speye(N^2)*eps^2/(dt^2*a22^2) - (2+d)/d*Lp);
    RHS = eps^2/(dt^2*a22^2)*p_2 - 1/(dt*a22)*(2+d)/d*Div_O2(u1_2,u2_2,dx,dy);
    % Solve
    rhs = reshape(RHS,[],1); 
    p_ = Lp2\rhs; 
    p_2 = reshape(p_,N,N); 

  
    %% Compute u^n+1
    comp11 = a22*Theta11_2;
    comp12 = a22*Theta12_2;
    comp21 = a22*Theta21_2;
    comp22 = a22*Theta22_2;
    flux_u1_2LF = Lax_Friedrichs3(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_2LF = Lax_Friedrichs3(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx2,p_dy2] = Grad_O2(p_2,dx,dy);

    u1_2 = u1 - dt*flux_u1_2LF  - dt*a22*p_dx2;
    u2_2 = u2 - dt*flux_u2_2LF  - dt*a22*p_dy2;


    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_2 = rho + d/(d+2)*(eps*p_2 - eps*p_2star);
    T_2 = T_2 + 2/d*d/(d+2)*(eps*p_2 - eps*p_2star);

    % %% Compute rho^n+1
    % u_dx_rho2 = Lax_Friedrichs(a22*u1_2/eps,a22*u2_2/eps,rho,alpha_e,dx,dy);
    % rho_2 = rho - dt*u_dx_rho2;%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute T^n+1
    % u_dx_T2 = Lax_Friedrichs(a22*u1_2/eps,a22*u2_2/eps,T,alpha_e,dx,dy);
    % T_2 = T_2 - dt*2/d*u_dx_T2;%Div_O2(u1,u2,dx,dy);


    %% --------------- STAGE 3 -------------------------------------------
    %% Compute Gamma^n+1
    [T_dx2,T_dy2] = Grad_O2(T_2,dx,dy);

    Gamma1_3 = eps^2*tau/(eps^2*tau+dt*a33)*Gamma1 + ...
        (d+2)*dt/(eps^2*tau+dt*a33)*(a31_*(u1_1.*T_1 - T_dx1*tau)...
        +a32_*(u1_2.*T_2 - T_dx2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Gamma1_2;
    Gamma2_3 = eps^2*tau/(eps^2*tau+dt*a33)*Gamma2 + ...
        (d+2)*dt/(eps^2*tau+dt*a33)*(a31_*(u2_1.*T_1 - T_dy1*tau)...
        +a32_*(u2_2.*T_2 - T_dy2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Gamma2_2;

    %% Compute p^n+1/2 and T^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_Gamma1_3 = Div_O2(Gamma1_3,Gamma2_3,dx,dy);%Lax_Friedrichs1(Gamma1_3,Gamma2_3,T,alpha_i,dx,dy);
    flux_Gamma2_3 = flux_Gamma1_3;% Lax_Friedrichs1(Gamma1_3,Gamma2_3,p,alpha_i,dx,dy);
    p_3 = p - dt/d*(a32*flux_Gamma2_2 + a33*flux_Gamma2_3)/eps;
    p_3star = p_3; % store value
    T_3 = T - dt/d*(a32*flux_Gamma1_2 + a33*flux_Gamma1_3);


    %% Compute Theta^n+1
    Du11_2 = Div_O2(2*u1_2,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_2,0*u1,Theta11,alpha_i,dx,dy);
    Du12_2 = Div_O2(u2_2,u1_2,dx,dy);%Lax_Friedrichs1(u2_2,u1_2,Theta12,alpha_i,dx,dy);
    Du21_2 = Div_O2(u2_2,u1_2,dx,dy);%Lax_Friedrichs1(u2_2,u1_2,Theta21,alpha_i,dx,dy);
    Du22_2 = Div_O2(u2*0,2*u2_2,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_2,Theta22,alpha_i,dx,dy);


    Au11_2 = u1_2.^2 - (u1_2.^2+u2_2.^2)/2;
    Au12_2 = u1_2.*u2_2;
    Au21_2 = u1_2.*u2_2;
    Au22_2 = u2_2.^2 - (u1_2.^2+u2_2.^2)/2;


    div_u_2 = Div_O2(u1_2,u2_2,dx,dy);


    Theta11_3 = eps^2*tau/(eps^2*tau+dt*a33)*Theta11 +...
        +dt/(eps^2*tau+dt*a33)*(a31_*(Au11_1 - Du11_1*tau)...
        +a32_*(Au11_2 - Du11_2*tau)) ...
        + dt*tau/(eps^2*tau+dt*a33)*2/d*(a31_*div_u_1 + a32_*div_u_2)...
        - dt/(eps^2*tau+dt*a33)*a32*Theta11_2;

    Theta12_3 = eps^2*tau/(eps^2*tau+dt*a33)*Theta12 +...
        +dt/(eps^2*tau+dt*a33)*(a31_*(Au12_1 - Du12_1*tau)...
        + a32_*(Au12_2 - Du12_2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Theta12_2;

    Theta21_3 = eps^2*tau/(eps^2*tau+dt*a33)*Theta21 +...
        +dt/(eps^2*tau+dt*a33)*(a31_*(Au21_1 - Du21_1*tau)...
        + a32_*(Au21_2 - Du21_2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Theta21_2;

    Theta22_3 = eps^2*tau/(eps^2*tau+dt*a33)*Theta22 +...
        +dt/(eps^2*tau+dt*a33)*(a31_*(Au22_1 - Du22_1*tau)...
        + a32_*(Au22_2 - Du22_2*tau)) ...
        + dt*tau/(eps^2*tau+dt*a33)*2/d*(a31_*div_u_1 + a32_*div_u_2)...
        - dt/(eps^2*tau+dt*a33)*a32*Theta22_2;


    %% p-Form
    %% Compute u^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_u1_3 = Div_O2(Theta11_3,Theta12_3,dx,dy);%Lax_Friedrichs1(Theta11_3,Theta12_3,u1,alpha_i,dx,dy);
    flux_u2_3 = Div_O2(Theta21_3,Theta22_3,dx,dy);%Lax_Friedrichs1(Theta21_3,Theta22_3,u2,alpha_i,dx,dy);

    u1_3 = u1 - dt*(a32*flux_u1_2 + a33*flux_u1_3);
    u2_3 = u2 - dt*(a32*flux_u2_2 + a33*flux_u2_3);


    %% Compute p^n+1   
    Lp3 = (speye(N^2)*eps^2/(dt^2*a33^2) - (2+d)/d*Lp);

    D_u2 = Div_O2(u1_2,u2_2,dx,dy); 
    Lapl_p2 = Lapl_O2(p_2,dx,dy);
    RHS = eps^2/(dt^2*a33^2)*p_3 - 1/(dt*a33)*(2+d)/d*Div_O2(u1_3,u2_3,dx,dy)...
        - 1/(dt*a33^2)*a32*(2+d)/d*D_u2 + 1/a33*a32*(2+d)/d*Lapl_p2;
    % Solve
    rhs = reshape(RHS,[],1); 
    p_ = Lp3\rhs; 

    p_3 = reshape(p_,N,N); 

  
    %% Compute u^n+1
    comp11 = a32*Theta11_2 + a33*Theta11_3;
    comp12 = a32*Theta12_2 + a33*Theta12_3;
    comp21 = a32*Theta21_2 + a33*Theta21_3;
    comp22 = a32*Theta22_2 + a33*Theta22_3;
    flux_u1_3LF = Lax_Friedrichs3(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_3LF = Lax_Friedrichs3(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx3,p_dy3] = Grad_O2(p_3,dx,dy);

    u1_3 = u1 - dt*flux_u1_3LF  - dt*(a32*p_dx2 + a33*p_dx3);
    u2_3 = u2 - dt*flux_u2_3LF   - dt*(a32*p_dy2 + a33*p_dy3);

    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_3 = rho + d/(d+2)*(eps*p_3 - eps*p_3star);
    T_3 = T_3 + 2/d*d/(d+2)*(eps*p_3 - eps*p_3star);


    % %% Compute rho^n+1
    % u_dx_rho3 = Lax_Friedrichs((a32*u1_2+a33*u1_3),(a32*u2_2+a33*u2_3),...
    %     rho,alpha_e,dx,dy);
    % rho_3 = rho - dt*(u_dx_rho3);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute T^n+1
    % u_dx_T3 = Lax_Friedrichs((a32*u1_2+a33*u1_3),(a32*u2_2+a33*u2_3),...
    %     T,alpha_e,dx,dy);
    % T_3 = T_3 - dt*2/d*(u_dx_T3);%Div_O2(u1,u2,dx,dy);

    %% UPDATE VARS
    T = T_3; rho = rho_3; u1 = u1_3; u2 = u2_3; p = p_3; 
    Gamma1 = Gamma1_3; Gamma2 = Gamma2_3;
    Theta11 = Theta11_3; Theta12 = Theta12_3; Theta21 = Theta21_3;
    Theta22 = Theta22_3;
    %% ________PLOT_______________________________
    u2_h(:,2:n+1) = u2;
    u2_h(:,1) = u2(:,end-1);
    u2_h(:,n+2) = u2(:,2);

    u1_h(2:n+1,:) = u1;
    u1_h(1,:) = u1(end-1,:);
    u1_h(n+2,:) = u1(2,:);
    du1_dy = 1/(2*dx)*(u1_h(3:end,:) - u1_h(1:end-2,:));
    du2_dx = 1/(2*dx)*(u2_h(:,3:end) - u2_h(:,1:end-2));
    vort = du2_dx - du1_dy;

    contour(X,Y,vort,30)
    title(['t is ',num2str(t)])
    drawnow
end
toc;
%% SAVE DATA

    % savename = ['Data_2D_ConvTest_3rdPaper/RK2WENO3O2LFalgFin','_tau=0','eps=1e-6','_',int2str(n),'.mat'];
    % % savename = ['1DTests-IMEX_RK2/RK2kWENO3O2_tau=01','eps=1','M=',int2str(M),'n=',int2str(n)];
    % save(savename, 'X','Y','dx','rho','u1','u2','T','p');
% end



