%% Relaxationsys in 2D - Double Shear Layer test case
% Feb/26
% Relaxationssystem mit p Formulierung
% IMEX RK3 

% clear; close all; clc;
function Algo2DRelSysRK3(n)
eps = 1e-06; % bis 1e-01 stabil f. alle tau

tau = 0;
delta = 0.05;
rho = pi/15;
% Skaliere die Schichtdicke: pi/15 im [0,2pi]-Gebiet entspricht 1/30 im [0,1]-Gebiet
% rho = 1/30;
tFinal = 1;
% tFinal = 9;
CFL = 0.75;

d = 2; % dimension

% n = 150; % number of grid pnts


% IMEX Verfahren - Typ A - erste Zeile d. Tableaus 0

a21_ = 1/2; 
a31_ = 11/18; a32_ = 1/18;
a41_ = 5/6; a42_ = -5/6; a43_ = 1/2;
a51_ = 1/4; a52_ = 7/4; a53_ = 3/4; a54_ = -7/4;

a11 = 0;
a21 = 0; a22 = 1/2;
a31 = 0; a32 = 1/6; a33 = 1/2;
a41 = 0; a42 = -1/2; a43 = 1/2; a44 = 1/2;
a51 = 0; a52 = 3/2; a53 = -3/2; a54 = 1/2; a55 = 1/2;

% Discretize spatial domain
x_min = 0; x_max = 2*pi; 
x = linspace(x_min,x_max,n); y = x; dx = x(2)-x(1); dy = dx;

% meshgrid
[X,Y] = meshgrid(x,y); 

%% initial conditions
% initial condition like in Lattice-Boltzmann paper
u1 = (Y <= pi).*tanh(1/rho*(Y-pi/2)) + (Y> pi).*tanh(1/rho*(3*pi/2-Y));
% Skalierte Initial Condition: 
% pi -> 0.5 (Mitte)
% pi/2 -> 0.25 (untere Scherschicht)
% 3*pi/2 -> 0.75 (obere Scherschicht)
% u1 = (Y <= 0.5).*tanh(1/rho*(Y-0.25)) + (Y > 0.5).*tanh(1/rho*(0.75-Y));
u2 = delta*sin(X);
% Störung u2: Sinuswelle muss über das Gebiet [0,1] eine volle Periode haben (Daher * 2*pi)
% u2 = delta*sin(2*pi*X);


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
    r2 = zeros(1,N); r2(1) = -30; r2(2) = 16; r2(3) = -1; A2_ = toeplitz(r2); 
    A2_(1,N-1) = 16; A2_(N,2) = 16; A2_(2,N-1) = -1; A2_(N,3) = -1; A2_(1,N-2) = -1; A2_(N-1,2) =-1;
    Lp = (kron(speye(N),sparse(A2_)) + kron(sparse(A2_),speye(N)))/(12*dx^2);
    % Lp = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    alpha_e = 1; alpha_i = 0;
    opts.type = 'ilutp'; 
    opts.droptol = 1e-8; 

    %% --------------- STAGE 1 -------------------------------------------
    u1_1 = u1; u2_1 = u2; rho_1 = rho; p_1 = p; T_1 = T;
    Theta11_1 = Theta11; Theta12_1 = Theta12; Theta21_1 =  Theta21; Theta22_1 = Theta22;
    Gamma1_1 = Gamma1; Gamma2_1 = Gamma2;

    %% --------------- STAGE 2 -------------------------------------------
    %% Compute Gamma^n+1
    [T_dx1,T_dy1] = Grad_O4(T_1,dx,dy);

    Gamma1_2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma1 + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u1_1.*T_1 - T_dx1*tau);
    Gamma2_2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma2 + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u2_1.*T_1 - T_dy1*tau);

    %% Compute p^n+1/2 and T^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_Gamma1_2 = Div_O4(Gamma1_2,Gamma2_2,dx,dy);%Lax_Friedrichs1(Gamma1_2,Gamma2_2,T,alpha_i,dx,dy);
    flux_Gamma2_2 = flux_Gamma1_2;%Lax_Friedrichs1(Gamma1_2,Gamma2_2,p,alpha_i,dx,dy);
    p_2 = p - dt/d*a22*flux_Gamma2_2/eps;
    p_2star = p_2; % introduce helper var.
    T_2 = T - dt/d*a22*flux_Gamma1_2;


    %% Compute Theta^n+1
    Du11_1 = Div_O4(2*u1_1,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_1,0*u1,Theta11,alpha_i,dx,dy);
    Du12_1 = Div_O4(u2_1,u1_1,dx,dy);%Lax_Friedrichs1(u2_1,u1_1,Theta12,alpha_i,dx,dy);
    Du21_1 = Div_O4(u2_1,u1_1,dx,dy);%Lax_Friedrichs1(u2_1,u1_1,Theta21,alpha_i,dx,dy);
    Du22_1 = Div_O4(u2*0,2*u2_1,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_1,Theta22,alpha_i,dx,dy);


    Au11_1 = u1_1.^2 - (u1_1.^2+u2_1.^2)/2;
    Au12_1 = u1_1.*u2_1;
    Au21_1 = u1_1.*u2_1;
    Au22_1 = u2_1.^2 - (u1_1.^2+u2_1.^2)/2;


    div_u_1 = Div_O4(u1_1,u2_1,dx,dy);


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
    flux_u1_2 = Div_O4(Theta11_2,Theta12_2,dx,dy);%Lax_Friedrichs1(Theta11_2,Theta12_2,u1,alpha_i,dx,dy);
    flux_u2_2 = Div_O4(Theta21_2,Theta22_2,dx,dy);%Lax_Friedrichs1(Theta21_2,Theta22_2,u2,alpha_i,dx,dy);

    u1_2 = u1 - dt*a22*flux_u1_2;
    u2_2 = u2 - dt*a22*flux_u2_2;


    %% Compute p^n+1   
    Lp2 = (speye(N^2)*eps^2/(dt^2*a22^2) - (2+d)/d*Lp);
    RHS = eps^2/(dt^2*a22^2)*p_2 - 1/(dt*a22)*(2+d)/d*Div_O4(u1_2,u2_2,dx,dy);
    % Solve
    rhs = reshape(RHS,[],1); 
    
    % [L, U] = ilu(Lp2, opts);
    p_ = Lp2\rhs; 
    % p_ = gmres(Lp2, rhs, 100,1e-12, 200, L, U);
    p_2 = reshape(p_,N,N); 


    %% Compute u^n+1
    comp11 = a22*Theta11_2;
    comp12 = a22*Theta12_2;
    comp21 = a22*Theta21_2;
    comp22 = a22*Theta22_2;
    flux_u1_2LF = Lax_Friedrichs(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_2LF = Lax_Friedrichs(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx2,p_dy2] = Grad_O4(p_2,dx,dy);
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
    [T_dx2,T_dy2] = Grad_O4(T_2,dx,dy);

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
    flux_Gamma1_3 = Div_O4(Gamma1_3,Gamma2_3,dx,dy);%Lax_Friedrichs1(Gamma1_3,Gamma2_3,T,alpha_i,dx,dy);
    flux_Gamma2_3 = flux_Gamma1_3;%Lax_Friedrichs1(Gamma1_3,Gamma2_3,p,alpha_i,dx,dy);
    p_3 = p - dt/d*(a32*flux_Gamma2_2 + a33*flux_Gamma2_3)/eps;
    p_3star = p_3;
    T_3 = T - dt/d*(a32*flux_Gamma1_2 + a33*flux_Gamma1_3);


    %% Compute Theta^n+1
    Du11_2 = Div_O4(2*u1_2,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_2,0*u1,Theta11,alpha_i,dx,dy);
    Du12_2 = Div_O4(u2_2,u1_2,dx,dy);%Lax_Friedrichs1(u2_2,u1_2,Theta12,alpha_i,dx,dy);
    Du21_2 = Div_O4(u2_2,u1_2,dx,dy);%Lax_Friedrichs1(u2_2,u1_2,Theta21,alpha_i,dx,dy);
    Du22_2 = Div_O4(u2*0,2*u2_2,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_2,Theta22,alpha_i,dx,dy);


    Au11_2 = u1_2.^2 - (u1_2.^2+u2_2.^2)/2;
    Au12_2 = u1_2.*u2_2;
    Au21_2 = u1_2.*u2_2;
    Au22_2 = u2_2.^2 - (u1_2.^2+u2_2.^2)/2;


    div_u_2 = Div_O4(u1_2,u2_2,dx,dy);


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
    flux_u1_3 = Div_O4(Theta11_3,Theta12_3,dx,dy);%Lax_Friedrichs1(Theta11_3,Theta12_3,u1,alpha_i,dx,dy);
    flux_u2_3 = Div_O4(Theta21_3,Theta22_3,dx,dy);%Lax_Friedrichs1(Theta21_3,Theta22_3,u2,alpha_i,dx,dy);

    u1_3 = u1 - dt*(a32*flux_u1_2 + a33*flux_u1_3);
    u2_3 = u2 - dt*(a32*flux_u2_2 + a33*flux_u2_3);


    %% Compute p^n+1   
    Lp3 = (speye(N^2)*eps^2/(dt^2*a33^2) - (2+d)/d*Lp);
    D_u2 = Div_O4(u1_2,u2_2,dx,dy); Lapl_p2 = Lapl_O4(p_2,dx,dy);
    RHS = eps^2/(dt^2*a33^2)*p_3 - 1/(dt*a33)*(2+d)/d*Div_O4(u1_3,u2_3,dx,dy)...
        - 1/(dt*a33^2)*a32*(2+d)/d*D_u2 + 1/a33*a32*(2+d)/d*Lapl_p2;
    % Solve
    rhs = reshape(RHS,[],1); 
    p_ = Lp3\rhs; 
    % [L, U] = ilu(Lp3, opts);
    % p_ = gmres(Lp3, rhs, 100,1e-12, 200, L, U);
    p_3 = reshape(p_,N,N); 


    %% Compute u^n+1
    comp11 = a32*Theta11_2 + a33*Theta11_3;
    comp12 = a32*Theta12_2 + a33*Theta12_3;
    comp21 = a32*Theta21_2 + a33*Theta21_3;
    comp22 = a32*Theta22_2 + a33*Theta22_3;
    flux_u1_3LF = Lax_Friedrichs(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_3LF = Lax_Friedrichs(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx3,p_dy3] = Grad_O4(p_3,dx,dy);

    u1_3 = u1 - dt*flux_u1_3LF  - dt*(a32*p_dx2 + a33*p_dx3);
    u2_3 = u2 - dt*flux_u2_3LF   - dt*(a32*p_dy2 + a33*p_dy3);

    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_3 = rho + d/(d+2)*(eps*p_3 - eps*p_3star);
    T_3 = T_3 + 2/d*d/(d+2)*(eps*p_3 - eps*p_3star);

    % %% Compute rho^n+1
    % u_dx_rho3 = Lax_Friedrichs((a32*u1_2+a33*u1_3)/eps,(a32*u2_2+a33*u2_3)/eps,...
    %     rho,alpha_e,dx,dy);
    % rho_3 = rho - dt*(u_dx_rho3);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute T^n+1
    % u_dx_T3 = Lax_Friedrichs((a32*u1_2+a33*u1_3)/eps,(a32*u2_2+a33*u2_3)/eps,...
    %     T,alpha_e,dx,dy);
    % T_3 = T_3 - dt*2/d*(u_dx_T3);%Div_O2(u1,u2,dx,dy);


    %% --------------- STAGE 4 -------------------------------------------
    %% Compute Gamma^n+1
    [T_dx3,T_dy3] = Grad_O4(T_3,dx,dy);

    Gamma1_4 = eps^2*tau/(eps^2*tau+dt*a44)*Gamma1 + ...
        (d+2)*dt/(eps^2*tau+dt*a44)*(a41_*(u1_1.*T_1 - T_dx1*tau)...
        +a42_*(u1_2.*T_2 - T_dx2*tau) +a43_*(u1_3.*T_3 - T_dx3*tau))...
        - dt/(eps^2*tau+dt*a44)*(a42*Gamma1_2 + a43*Gamma1_3);
    Gamma2_4 = eps^2*tau/(eps^2*tau+dt*a44)*Gamma2 + ...
        (d+2)*dt/(eps^2*tau+dt*a44)*(a41_*(u2_1.*T_1 - T_dy1*tau)...
        +a42_*(u2_2.*T_2 - T_dy2*tau) +a43_*(u2_3.*T_3 - T_dy3*tau))...
        - dt/(eps^2*tau+dt*a44)*(a42*Gamma2_2 + a43*Gamma2_3);

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma1_4 = Div_O4(Gamma1_4,Gamma2_4,dx,dy);%Lax_Friedrichs1(Gamma1_4,Gamma2_4,T,alpha_i,dx,dy);
    flux_Gamma2_4 = flux_Gamma1_4;%Lax_Friedrichs1(Gamma1_4,Gamma2_4,p,alpha_i,dx,dy);
    p_4 = p - dt/d*(a42*flux_Gamma2_2 + a43*flux_Gamma2_3 + a44*flux_Gamma2_4)/eps;
    p_4star = p_4;
    T_4 = T - dt/d*(a42*flux_Gamma1_2 + a43*flux_Gamma1_3 + a44*flux_Gamma1_4);


    %% Compute Theta^n+1
    Du11_3 = Div_O4(2*u1_3,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_3,0*u1,Theta11,alpha_i,dx,dy);
    Du12_3 = Div_O4(u2_3,u1_3,dx,dy);%Lax_Friedrichs1(u2_3,u1_3,Theta12,alpha_i,dx,dy);
    Du21_3 = Div_O4(u2_3,u1_3,dx,dy);%Lax_Friedrichs1(u2_3,u1_3,Theta21,alpha_i,dx,dy);
    Du22_3 = Div_O4(u2*0,2*u2_3,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_3,Theta22,alpha_i,dx,dy);


    Au11_3 = u1_3.^2 - (u1_3.^2+u2_3.^2)/2;
    Au12_3 = u1_3.*u2_3;
    Au21_3 = u1_3.*u2_3;
    Au22_3 = u2_3.^2 - (u1_3.^2+u2_3.^2)/2;


    div_u_3 = Div_O4(u1_3,u2_3,dx,dy);


    Theta11_4 = eps^2*tau/(eps^2*tau+dt*a44)*Theta11 +...
        +dt/(eps^2*tau+dt*a44)*(a41_*(Au11_1 - Du11_1*tau)...
        +a42_*(Au11_2 - Du11_2*tau) +a43_*(Au11_3 - Du11_3*tau)) ...
        + dt*tau/(eps^2*tau+dt*a44)*2/d*(a41_*div_u_1 + a42_*div_u_2 + a43_*div_u_3)...
        - dt/(eps^2*tau+dt*a44)*(a42*Theta11_2 + a43*Theta11_3);

    Theta12_4 = eps^2*tau/(eps^2*tau+dt*a44)*Theta12 +...
        +dt/(eps^2*tau+dt*a44)*(a41_*(Au12_1 - Du12_1*tau)...
        + a42_*(Au12_2 - Du12_2*tau) + a43_*(Au12_3 - Du12_3*tau))...
        - dt/(eps^2*tau+dt*a44)*(a42*Theta12_2 + a43*Theta12_3);

    Theta21_4 = eps^2*tau/(eps^2*tau+dt*a44)*Theta21 +...
        +dt/(eps^2*tau+dt*a44)*(a41_*(Au21_1 - Du21_1*tau)...
        + a42_*(Au21_2 - Du21_2*tau) + a43_*(Au21_3 - Du21_3*tau))...
        - dt/(eps^2*tau+dt*a44)*(a42*Theta21_2 +a43*Theta21_3);

    Theta22_4 = eps^2*tau/(eps^2*tau+dt*a44)*Theta22 +...
        +dt/(eps^2*tau+dt*a44)*(a41_*(Au22_1 - Du22_1*tau)...
        + a42_*(Au22_2 - Du22_2*tau) + a43_*(Au22_3 - Du22_3*tau)) ...
        + dt*tau/(eps^2*tau+dt*a44)*2/d*(a41_*div_u_1 + a42_*div_u_2 + a43_*div_u_3)...
        - dt/(eps^2*tau+dt*a44)*(a42*Theta22_2 +a43*Theta22_3);


    %% p-Form
    %% Compute u^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_u1_4 = Div_O4(Theta11_4,Theta12_4,dx,dy);%Lax_Friedrichs1(Theta11_4,Theta12_4,u1,alpha_i,dx,dy);
    flux_u2_4 = Div_O4(Theta21_4,Theta22_4,dx,dy);%Lax_Friedrichs1(Theta21_4,Theta22_4,u2,alpha_i,dx,dy);

    u1_4 = u1 - dt*(a42*flux_u1_2 + a43*flux_u1_3 + a44*flux_u1_4);
    u2_4 = u2 - dt*(a42*flux_u2_2 + a43*flux_u2_3 + a44*flux_u2_4);


    %% Compute p^n+1   
    Lp4 = (speye(N^2)*eps^2/(dt^2*a44^2) - (2+d)/d*Lp);
    % D_u2 = Div_O2(u1_2,u2_2,dx,dy); Lapl_p2 = Lapl_O2(p_2,dx,dy); 
    Lapl_p3 = Lapl_O4(p_3,dx,dy); D_u3 = Div_O4(u1_3,u2_3,dx,dy);
    RHS = eps^2/(dt^2*a44^2)*p_4 - 1/(dt*a44)*(2+d)/d*Div_O4(u1_4,u2_4,dx,dy)...
        - 1/(dt*a44^2)*a42*(2+d)/d*D_u2 - 1/(dt*a44^2)*a43*(2+d)/d*D_u3 ...
        + 1/a44*a42*(2+d)/d*Lapl_p2 + 1/a44*a43*(2+d)/d*Lapl_p3;
    % Solve
    rhs = reshape(RHS,[],1); 
    p_ = Lp4\rhs; 
    % [L, U] = ilu(Lp4, opts);
    % p_ = gmres(Lp4, rhs, 100,1e-12, 200, L, U);
    p_4 = reshape(p_,N,N); 


    %% Compute u^n+1
    comp11 = a42*Theta11_2 + a43*Theta11_3 + a44*Theta11_4;
    comp12 = a42*Theta12_2 + a43*Theta12_3 + a44*Theta12_4;
    comp21 = a42*Theta21_2 + a43*Theta21_3 + a44*Theta21_4;
    comp22 = a42*Theta22_2 + a43*Theta22_3 + a44*Theta22_4;
    flux_u1_4LF = Lax_Friedrichs(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_4LF = Lax_Friedrichs(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx4,p_dy4] = Grad_O4(p_4,dx,dy);
    u1_4 = u1 - dt*flux_u1_4LF  - dt*(a42*p_dx2 + a43*p_dx3 + a44*p_dx4);
    u2_4 = u2 - dt*flux_u2_4LF  - dt*(a42*p_dy2 + a43*p_dy3 + a44*p_dy4);

    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_4 = rho + d/(d+2)*(eps*p_4 - eps*p_4star);
    T_4 = T_4 + 2/d*d/(d+2)*(eps*p_4 - eps*p_4star);

    % %% Compute rho^n+1
    % u_dx_rho4 = Lax_Friedrichs((a42*u1_2+a43*u1_3+a44*u1_4)/eps,...
    %     (a42*u2_2+a43*u2_3+a44*u2_4)/eps,...
    %     rho,alpha_e,dx,dy);
    % rho_4 = rho - dt*(u_dx_rho4);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute T^n+1
    % u_dx_T4 = Lax_Friedrichs((a42*u1_2+a43*u1_3+a44*u1_4)/eps,...
    %     (a42*u2_2+a43*u2_3+a44*u2_4)/eps,...
    %     T,alpha_e,dx,dy);
    % T_4 = T_4 - dt*2/d*(u_dx_T4);%Div_O2(u1,u2,dx,dy);


    %% --------------- STAGE 5 -------------------------------------------
    %% Compute Gamma^n+1
    [T_dx4,T_dy4] = Grad_O4(T_4,dx,dy);

    Gamma1_5 = eps^2*tau/(eps^2*tau+dt*a55)*Gamma1 + ...
        (d+2)*dt/(eps^2*tau+dt*a55)*(a51_*(u1_1.*T_1 - T_dx1*tau)...
        +a52_*(u1_2.*T_2 - T_dx2*tau) +a53_*(u1_3.*T_3 - T_dx3*tau)...
        +a54_*(u1_4.*T_4 - T_dx4*tau))...
        - dt/(eps^2*tau+dt*a55)*(a52*Gamma1_2 + a53*Gamma1_3 + a54*Gamma1_4);
    Gamma2_5 = eps^2*tau/(eps^2*tau+dt*a55)*Gamma2 + ...
        (d+2)*dt/(eps^2*tau+dt*a55)*(a51_*(u2_1.*T_1 - T_dy1*tau)...
        +a52_*(u2_2.*T_2 - T_dy2*tau) +a53_*(u2_3.*T_3 - T_dy3*tau)...
        +a54_*(u2_4.*T_4 - T_dy4*tau))...
        - dt/(eps^2*tau+dt*a55)*(a52*Gamma2_2 + a53*Gamma2_3 + a54*Gamma2_4);

    %% Compute p^n+1/2 and T^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_Gamma1_5 = Div_O4(Gamma1_5,Gamma2_5,dx,dy);%Lax_Friedrichs1(Gamma1_4,Gamma2_4,T,alpha_i,dx,dy);
    flux_Gamma2_5 = flux_Gamma1_5;%Lax_Friedrichs1(Gamma1_4,Gamma2_4,p,alpha_i,dx,dy);
    p_5 = p - dt/d*(a52*flux_Gamma2_2 + a53*flux_Gamma2_3...
        + a54*flux_Gamma2_4 +a55*flux_Gamma2_5)/eps;
    p_5star = p_5;
    T_5 = T - dt/d*(a52*flux_Gamma1_2 + a53*flux_Gamma1_3 ...
        + a54*flux_Gamma1_4 + a55*flux_Gamma1_5);


    %% Compute Theta^n+1
    Du11_4 = Div_O4(2*u1_4,0*u1,dx,dy);%Lax_Friedrichs1(2*u1_4,0*u1,Theta11,alpha_i,dx,dy);
    Du12_4 = Div_O4(u2_4,u1_4,dx,dy);%Lax_Friedrichs1(u2_4,u1_4,Theta12,alpha_i,dx,dy);
    Du21_4 = Div_O4(u2_4,u1_4,dx,dy);%Lax_Friedrichs1(u2_4,u1_4,Theta21,alpha_i,dx,dy);
    Du22_4 = Div_O4(u2*0,2*u2_4,dx,dy);%Lax_Friedrichs1(u2*0,2*u2_4,Theta22,alpha_i,dx,dy);


    Au11_4 = u1_4.^2 - (u1_4.^2+u2_4.^2)/2;
    Au12_4 = u1_4.*u2_4;
    Au21_4 = u1_4.*u2_4;
    Au22_4 = u2_4.^2 - (u1_4.^2+u2_4.^2)/2;


    div_u_4 = Div_O4(u1_4,u2_4,dx,dy);


    Theta11_5 = eps^2*tau/(eps^2*tau+dt*a55)*Theta11 +...
        +dt/(eps^2*tau+dt*a55)*(a51_*(Au11_1 - Du11_1*tau)...
        +a52_*(Au11_2 - Du11_2*tau) +a53_*(Au11_3 - Du11_3*tau)...
        +a54_*(Au11_4 - Du11_4*tau)) ...
        + dt*tau/(eps^2*tau+dt*a55)*2/d*...
        (a51_*div_u_1 + a52_*div_u_2 + a53_*div_u_3 + a54_*div_u_4)...
        - dt/(eps^2*tau+dt*a55)*(a52*Theta11_2 + a53*Theta11_3 + a54*Theta11_4);

    Theta12_5 = eps^2*tau/(eps^2*tau+dt*a55)*Theta12 +...
        +dt/(eps^2*tau+dt*a55)*(a51_*(Au12_1 - Du12_1*tau)...
        + a52_*(Au12_2 - Du12_2*tau) + a53_*(Au12_3 - Du12_3*tau)...
         + a54_*(Au12_4 - Du12_4*tau))...
        - dt/(eps^2*tau+dt*a55)*(a52*Theta12_2 + a53*Theta12_3 + a54*Theta12_4);

    Theta21_5 = eps^2*tau/(eps^2*tau+dt*a55)*Theta21 +...
        +dt/(eps^2*tau+dt*a55)*(a51_*(Au21_1 - Du21_1*tau)...
        + a52_*(Au21_2 - Du21_2*tau) + a53_*(Au21_3 - Du21_3*tau)...
         + a54_*(Au21_4 - Du21_4*tau))...
        - dt/(eps^2*tau+dt*a55)*(a52*Theta21_2 +a53*Theta21_3 +a54*Theta21_4);

    Theta22_5 = eps^2*tau/(eps^2*tau+dt*a55)*Theta22 +...
        +dt/(eps^2*tau+dt*a55)*(a51_*(Au22_1 - Du22_1*tau)...
        + a52_*(Au22_2 - Du22_2*tau) + a53_*(Au22_3 - Du22_3*tau)...
        + a54_*(Au22_4 - Du22_4*tau)) ...
        + dt*tau/(eps^2*tau+dt*a55)*2/d*...
        (a51_*div_u_1 + a52_*div_u_2 + a53_*div_u_3 + a54_*div_u_4)...
        - dt/(eps^2*tau+dt*a55)*(a52*Theta22_2 +a53*Theta22_3 +a54*Theta22_4);


    %% p-Form
    %% Compute u^n+1/2
    % For large eps, use alpha_e in the sollowing 2 lines
    flux_u1_5 = Div_O4(Theta11_5,Theta12_5,dx,dy);%Lax_Friedrichs1(Theta11_5,Theta12_5,u1,alpha_i,dx,dy);
    flux_u2_5 = Div_O4(Theta21_5,Theta22_5,dx,dy);%Lax_Friedrichs1(Theta21_5,Theta22_5,u2,alpha_i,dx,dy);

    u1_5 = u1 - dt*(a52*flux_u1_2 + a53*flux_u1_3 + a54*flux_u1_4 + a55*flux_u1_5);
    u2_5 = u2 - dt*(a52*flux_u2_2 + a53*flux_u2_3 + a54*flux_u2_4 + a55*flux_u2_5);


    %% Compute p^n+1   
    Lp5 = (speye(N^2)*eps^2/(dt^2*a55^2) - (2+d)/d*Lp);
    % D_u2 = Div_O2(u1_2,u2_2,dx,dy); Lapl_p2 = Lapl_O2(p_2,dx,dy); 
    % Lapl_p3 = Lapl_O2(p_3,dx,dy); D_u3 = Div_O2(u1_3,u2_3,dx,dy);
    Lapl_p4 = Lapl_O4(p_4,dx,dy); D_u4 = Div_O4(u1_4,u2_4,dx,dy);
    RHS = eps^2/(dt^2*a55^2)*p_5 - 1/(dt*a55)*(2+d)/d*Div_O4(u1_5,u2_5,dx,dy)...
        - 1/(dt*a55^2)*a52*(2+d)/d*D_u2 - 1/(dt*a55^2)*a53*(2+d)/d*D_u3 ...
        - 1/(dt*a55^2)*a54*(2+d)/d*D_u4 ...
        + 1/a55*a52*(2+d)/d*Lapl_p2 + 1/a55*a53*(2+d)/d*Lapl_p3 + 1/a55*a54*(2+d)/d*Lapl_p4;
    % Solve
    rhs = reshape(RHS,[],1); 
    p_ = Lp5\rhs; 
    % [L, U] = ilu(Lp5, opts);
    % p_ = gmres(Lp5, rhs, 100,1e-12, 200, L, U);
    p_5 = reshape(p_,N,N); 


    %% Compute u^n+1
    comp11 = a52*Theta11_2 + a53*Theta11_3 + a54*Theta11_4 + a55*Theta11_5;
    comp12 = a52*Theta12_2 + a53*Theta12_3 + a54*Theta12_4 + a55*Theta12_5;
    comp21 = a52*Theta21_2 + a53*Theta21_3 + a54*Theta21_4 + a55*Theta21_5;
    comp22 = a52*Theta22_2 + a53*Theta22_3 + a54*Theta22_4 + a55*Theta22_5;
    flux_u1_5LF = Lax_Friedrichs(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_5LF = Lax_Friedrichs(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx5,p_dy5] = Grad_O4(p_5,dx,dy);
    u1_5 = u1 - dt*flux_u1_5LF  - dt*(a52*p_dx2 + a53*p_dx3 + a54*p_dx4 + a55*p_dx5);
    u2_5 = u2 - dt*flux_u2_5LF  - dt*(a52*p_dy2 + a53*p_dy3 + a54*p_dy4 + a55*p_dy5);


    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_5 = rho + d/(d+2)*(eps*p_5 - eps*p_5star);
    T_5 = T_5 + 2/d*d/(d+2)*(eps*p_5 - eps*p_5star);

    % %% Compute rho^n+1
    % u_dx_rho5 = Lax_Friedrichs((a52*u1_2+a53*u1_3+a54*u1_4+a55*u1_5)/eps,...
    %     (a52*u2_2+a53*u2_3+a54*u2_4+a55*u2_5)/eps,...
    %     rho,alpha_e,dx,dy);
    % rho_5 = rho - dt*(u_dx_rho5);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute T^n+1
    % u_dx_T5 = Lax_Friedrichs((a52*u1_2+a53*u1_3+a54*u1_4+a55*u1_5)/eps,...
    %     (a52*u2_2+a53*u2_3+a54*u2_4+a55*u2_5)/eps,...
    %     T,alpha_e,dx,dy);
    % T_5 = T_5 - dt*2/d*(u_dx_T5);%Div_O2(u1,u2,dx,dy);

    %% UPDATE VARS
    T = T_5; rho = rho_5; u1 = u1_5; u2 = u2_5; p = p_5; 
    Gamma1 = Gamma1_5; Gamma2 = Gamma2_5;
    Theta11 = Theta11_5; Theta12 = Theta12_5; Theta21 = Theta21_5;
    Theta22 = Theta22_5;
    %% ________PLOT_______________________________
    u2_(:,2:n+1) = u2;
    u2_(:,1) = u2(:,end-1);
    u2_(:,n+2) = u2(:,2);

    u1_(2:n+1,:) = u1;
    u1_(1,:) = u1(end-1,:);
    u1_(n+2,:) = u1(2,:);
    du1_dy = 1/(2*dx)*(u1_(3:end,:) - u1_(1:end-2,:));
    du2_dx = 1/(2*dx)*(u2_(:,3:end) - u2_(:,1:end-2));
    vort = du2_dx - du1_dy;

    contour(X,Y,vort,30)
    title(['t is ',num2str(t)])
    drawnow
end

%% SAVE DATA
    savename = ['Data_2D_ConvTest_3rdPaper/RK3WENO5O4LFalgFin','_tau=0','eps=1e-6','_',int2str(n),'.mat'];
    % savename = ['1DTests-IMEX_RK2/RK2kWENO3O2_tau=01','eps=1','M=',int2str(M),'n=',int2str(n)];
    save(savename, 'X','Y','dx','rho','u1','u2','T','p');
    % savename = ['Data_2D_ConvTest_3rdPaper/RK3WENO5O4LFalg_tau=0','eps=1e-6','_',int2str(n),];
    % % savename = ['1DTests-IMEX_RK2/RK2kWENO3O2_tau=01','eps=1','M=',int2str(M),'n=',int2str(n)];
    % save(savename, 'X','Y','dx','rho','u1','u2','T','p');
end