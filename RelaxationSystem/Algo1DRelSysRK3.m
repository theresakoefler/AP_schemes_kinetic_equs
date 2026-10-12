%% Relaxationsys in 1D - fin
% Feb/26
% Relaxationssystem mit p Formulierung
% IMEX RK3 -> LaxFriedrichs1D() öffnen und auf WENO5 ändern!

clear; close all; clc;
eps = 1e-0;

tau = 0.001;
tFinal = 0.05;
CFL = 0.75;

d = 1; % dimension

n = 200; % number of grid pnts


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
x_min = 0; x_max = 1; x = linspace(x_min,x_max,n); y = x; dx = x(2)-x(1); 

%% initial conditions
u = x*0;
Gamma = u;
% initialize rho, 
rho = ones(1,n); 
rho(0 <= x & x <= 0.2) = 1;
rho(0.2 < x & x <= 0.3) = 1+eps^2;
rho(0.3 < x & x <= 0.7) = 1;
rho(0.7 < x & x <= 0.8) = 1-eps^2;

% initialize T
T = u*0; p = rho+T; 

t = 0; dt = dx^2*CFL; 

% Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); A_(1,N-1) = 1; A_(N,2) = 1;
    r2 = zeros(1,N); r2(1) = -30; r2(2) = 16; r2(3) = -1; A2_ = toeplitz(r2); 
    A2_(1,N-1) = 16; A2_(N,2) = 16; A2_(2,N-1) = -1; A2_(N,3) = -1; A2_(1,N-2) = -1; A2_(N-1,2) =-1;
    A2_ = A2_/(12*dx^2);
    A_ = A_/dx^2;

while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    alpha_e = 1; alpha_i = 0;

    %% --------------- STAGE 1 -------------------------------------------
    u1 = u; rho1 = rho; p1 = p; T1 = T;
    Gamma1 = Gamma; 

    %% --------------- STAGE 2 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx1 = Dx_O4(T1,dx);

    Gamma2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u1.*T1 - T_dx1*tau);

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma2 = Dx_O4(Gamma2,dx);
    p2 = p - dt/d*a22*flux_Gamma2;
    T2 = T - dt/d*a22*flux_Gamma2;


    %% Compute p^n+1   
    Lp2 = (speye(N)*eps^2/(dt^2*a22^2) - (2+d)/d*A2_);
    RHS = eps^2/(dt^2*a22^2)*p2 - eps/(dt*a22)*(2+d)/d*Dx_O4(u,dx);
    % Solve
    rhs = RHS'; 
    p_ = Lp2\rhs; p2 = p_'; 


    %% Compute u^n+1
    p_dx2 = Lax_Friedrichs1D(p2,u,alpha_e,dx);
    u2 = u  - dt*a22*p_dx2/eps;

    %% Compute rho^n+1
    u_dx_rho2 = Lax_Friedrichs1D(u2,rho,alpha_e,dx);
    rho2 = rho - dt/eps*a22*u_dx_rho2;

    %% Compute T^n+1
    u_dx_T2 = Lax_Friedrichs1D(u2,T,alpha_e,dx);
    T2 = T2 - dt/eps*2/d*a22*u_dx_T2;


    %% --------------- STAGE 3 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx2 = Dx_O4(T2,dx);

    Gamma3 = eps^2*tau/(eps^2*tau+dt*a33)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a33)*(a31_*(u1.*T1 - T_dx1*tau)...
        +a32_*(u2.*T2 - T_dx2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Gamma2;

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma3 = Dx_O4(Gamma3,dx);
    p3 = p - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3);
    T3 = T - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3);


    %% Compute p^n+1   
    Lp3 = (speye(N)*eps^2/(dt^2*a33^2) - (2+d)/d*A2_);
    Dx_u2 = Dx_O4(u2,dx); Dxx_p2 = Dxx_O4(p2,dx);
    RHS = eps^2/(dt^2*a33^2)*p3 - eps/(dt*a33)*(2+d)/d*Dx_O4(u,dx)...
        - eps/(dt*a33^2)*a32*(2+d)/d*Dx_u2 + 1/a33*a32*(2+d)/d*Dxx_p2;
    % Solve
    rhs = RHS'; 
    p_ = Lp3\rhs; p3 = p_'; 


    %% Compute u^n+1
    p_dx3 = Lax_Friedrichs1D(p3,u,alpha_e,dx);
    u3 = u  - dt*(a32*p_dx2 + a33*p_dx3)/eps;

    %% Compute rho^n+1
    u_dx_rho3 = Lax_Friedrichs1D(u3,rho,alpha_e,dx);
    rho3 = rho - dt/eps*(a32*u_dx_rho2 + a33*u_dx_rho3);

    %% Compute T^n+1
    u_dx_T3 = Lax_Friedrichs1D(u3,T,alpha_e,dx);
    T3 = T3 - dt/eps*2/d*(a32*u_dx_T2 + a33*u_dx_T3);


    %% --------------- STAGE 4 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx3 = Dx_O4(T3,dx);

    Gamma4 = eps^2*tau/(eps^2*tau+dt*a44)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a44)*(a41_*(u1.*T1 - T_dx1*tau)...
        +a42_*(u2.*T2 - T_dx2*tau) +a43_*(u3.*T3 - T_dx3*tau))...
        - dt/(eps^2*tau+dt*a44)*(a42*Gamma2 + a43*Gamma3);

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma4 = Dx_O4(Gamma4,dx);
    p4 = p - dt/d*(a42*flux_Gamma2 + a43*flux_Gamma3 + a44*flux_Gamma4);
    T4 = T - dt/d*(a42*flux_Gamma2 + a43*flux_Gamma3 + a44*flux_Gamma4);


    %% Compute p^n+1   
    Lp4 = (speye(N)*eps^2/(dt^2*a44^2) - (2+d)/d*A2_);
    Dxx_p3 = Dxx_O4(p3,dx); Dx_u3 = Dx_O4(u3,dx);
    RHS = eps^2/(dt^2*a44^2)*p4 - eps/(dt*a44)*(2+d)/d*Dx_O4(u,dx)...
        - eps/(dt*a44^2)*a42*(2+d)/d*Dx_u2 - eps/(dt*a44^2)*a43*(2+d)/d*Dx_u3 ...
        + 1/a44*a42*(2+d)/d*Dxx_p2 + 1/a44*a43*(2+d)/d*Dxx_p3;
    % Solve
    rhs = RHS'; 
    p_ = Lp4\rhs; p4 = p_'; 


    %% Compute u^n+1
    p_dx4 = Lax_Friedrichs1D(p4,u,alpha_e,dx);
    u4 = u  - dt*(a42*p_dx2 + a43*p_dx3 + a44*p_dx4)/eps;

    %% Compute rho^n+1
    u_dx_rho4 = Lax_Friedrichs1D(u4,rho,alpha_e,dx);
    rho4 = rho - dt/eps*(a42*u_dx_rho2 + a43*u_dx_rho3 + a44*u_dx_rho4);

    %% Compute T^n+1
    u_dx_T4 = Lax_Friedrichs1D(u4,T,alpha_e,dx);
    T4 = T4 - dt/eps*2/d*(a42*u_dx_T2 + a43*u_dx_T3 + a44*u_dx_T4);


    %% --------------- STAGE 5 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx4 = Dx_O4(T4,dx);

    Gamma5 = eps^2*tau/(eps^2*tau+dt*a55)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a55)*(a51_*(u1.*T1 - T_dx1*tau)...
        +a52_*(u2.*T2 - T_dx2*tau) +a53_*(u3.*T3 - T_dx3*tau)...
        +a54_*(u4.*T4 - T_dx4*tau))...
        - dt/(eps^2*tau+dt*a55)*(a52*Gamma2 + a53*Gamma3 + a54*Gamma4);

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma5 = Dx_O4(Gamma4,dx);
    p5 = p - dt/d*(a52*flux_Gamma2 + a53*flux_Gamma3...
        + a54*flux_Gamma4 +a55*flux_Gamma5);
    T5 = T - dt/d*(a52*flux_Gamma2 + a53*flux_Gamma3 ...
        + a54*flux_Gamma4 + a55*flux_Gamma5);


    %% Compute p^n+1   
    Lp5 = (speye(N)*eps^2/(dt^2*a55^2) - (2+d)/d*A2_);
    Dxx_p4 = Dxx_O4(p4,dx); Dx_u4 = Dx_O4(u4,dx);
    RHS = eps^2/(dt^2*a55^2)*p5 - eps/(dt*a55)*(2+d)/d*Dx_O4(u,dx)...
        - eps/(dt*a55^2)*a52*(2+d)/d*Dx_u2 - eps/(dt*a55^2)*a53*(2+d)/d*Dx_u3 ...
        - eps/(dt*a55^2)*a54*(2+d)/d*Dx_u4 ...
        + 1/a55*a52*(2+d)/d*Dxx_p2 + 1/a55*a53*(2+d)/d*Dxx_p3 + 1/a55*a54*(2+d)/d*Dxx_p4;
    % Solve
    rhs = RHS'; 
    p_ = Lp5\rhs; p5 = p_'; 


    %% Compute u^n+1
    p_dx5 = Lax_Friedrichs1D(p5,u,alpha_e,dx);
    u5 = u  - dt*(a52*p_dx2 + a53*p_dx3 + a54*p_dx4 + a55*p_dx5)/eps;

    %% Compute rho^n+1
    u_dx_rho5 = Lax_Friedrichs1D(u5,rho,alpha_e,dx);
    rho5 = rho - dt/eps*(a52*u_dx_rho2 + a53*u_dx_rho3 + a54*u_dx_rho4 + a55*u_dx_rho5);

    %% Compute T^n+1
    u_dx_T5 = Lax_Friedrichs1D(u5,T,alpha_e,dx);
    T5 = T5 - dt/eps*2/d*(a52*u_dx_T2 + a53*u_dx_T3 + a54*u_dx_T4 + a55*u_dx_T5);

    %% UPDATE VARS
    T = T5; rho = rho5; u = u5; p = p5; 
    Gamma = Gamma5;
    %% ________PLOT_______________________________
    subplot(3,1,1)
    plot(x,rho);
    title('rho')

    subplot(3,1,2)
    plot(x,u);
    title('u')

    subplot(3,1,3)
    plot(x,T);
    title('T')


    sgtitle(['t is ',num2str(t)])
    drawnow
end