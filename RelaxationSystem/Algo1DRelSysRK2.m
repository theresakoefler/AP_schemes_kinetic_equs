%% Relaxationsys in 1D - fin
% Feb/26
% Relaxationssystem mit p Formulierung
% IMEX RK2

clear; close all; clc;
tic;
eps = 0.01; % Geht bis eps = 1e-0

tau = 0.1;
tFinal = 0.0017%0.01;
CFL = 0.75;

d = 1; % dimension

n = 1001; % number of grid pnts


% IMEX Verfahren - Typ A - erste Zeile d. Tableaus 0
c = 1-1/sqrt(2); gamma = (c-1/2)/(c-1);

a21_ = c; 
a31_ = 1-1/(2*c); a32_ = 1/(2*c);

a11 = 0;
a21 = 0; a22 = c;
a31 = 0; a32 = 1-gamma; a33 = gamma;

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
T = u*0; p = (rho+T)/eps; 

t = 0; dt = dx^2*CFL; 

% Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); A_(1,N-1) = 1; A_(N,2) = 1;
    A_ = A_/dx^2;

while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    alpha_e = 1; alpha_i = 0;

    %% --------------- STAGE 1 -------------------------------------------
    u1 = u; rho1 = rho; p1 = p; T1 = T; Gamma1 = Gamma;

    %% --------------- STAGE 2 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx1 = Dx(T1,dx);

    Gamma2 = eps^2*tau/(eps^2*tau+dt*a22)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a22)*a21_*(u1.*T1 - T_dx1*tau);

    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma2 = Dx(Gamma2,dx);
    % p2 = p - dt/d*a22*flux_Gamma2;
    % T2 = T - dt/d*a22*flux_Gamma2;
    T2 =T - dt*a22*flux_Gamma2;
    p2 = p - dt/d*a22*flux_Gamma2/eps;
    p2_star = p2;

    %% Compute p^n+1   
    Lp2 = (speye(N)*eps^2/(dt^2*a22^2) - (2+d)/d*A_);
    % RHS = eps^2/(dt^2*a22^2)*p2 - eps/(dt*a22)*(2+d)/d*Dx(u,dx);
    u_dx = Dx(u,dx);
    RHS = eps^2/(dt^2*a22^2)*p2 - 1/(dt*a22)*(d+2)/d*u_dx;

    % Solve
    rhs = RHS'; 
    p_ = Lp2\rhs; p2 = p_'; 

  
    %% Compute u^n+1
    % p_dx2 = Lax_Friedrichs1D(p2,u,alpha_e,dx);
    % u2 = u - dt*a22*p_dx2/eps;
    p_dx = Lax_Friedrichs1D(p2*a22,u,alpha_e,dx); % MIT LF!!!
    u2 = u - dt*p_dx;

    %% Compute rho^n+1
    u_dx_rho2 = Lax_Friedrichs1D(u2,rho,alpha_e,dx);
    rho2 = rho - dt/eps*a22*u_dx_rho2;
    % 
    % %% Compute T^n+1
    u_dx_T2 = Lax_Friedrichs1D(u2,T,alpha_e,dx);
    T2 = T2 - dt/eps*2/d*a22*u_dx_T2;

    % rho2 = rho + d/(d+2)*(p2 - p2_star)*eps;
    % T2 = T2 + 2/d*d/(d+2)*(p2 - p2_star)*eps;


    %% --------------- STAGE 3 -------------------------------------------
    %% Compute Gamma^n+1
    T_dx2 = Dx(T2,dx);

    Gamma3 = eps^2*tau/(eps^2*tau+dt*a33)*Gamma + ...
        (d+2)*dt/(eps^2*tau+dt*a33)*(a31_*(u1.*T1 - T_dx1*tau)...
        +a32_*(u2.*T2 - T_dx2*tau))...
        - dt/(eps^2*tau+dt*a33)*a32*Gamma2;


    %% Compute p^n+1/2 and T^n+1/2
    flux_Gamma3 = Dx(Gamma3,dx);
    % p3 = p - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3);
    % T3 = T - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3);
    T3 =T - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3);
    p3 = p - dt/d*(a32*flux_Gamma2 + a33*flux_Gamma3)/eps;
    p3_star = p3;


    %% Compute p^n+1   
    Lp3 = (speye(N)*eps^2/(dt^2*a33^2) - (2+d)/d*A_);
    % Dx_u2 = Dx(u2,dx); Dxx_p2 = Dxx(p2,dx);
    % RHS = eps^2/(dt^2*a33^2)*p3 - eps/(dt*a33)*(2+d)/d*Dx(u,dx)...
    %     - eps/(dt*a33^2)*a32*(2+d)/d*Dx_u2 + 1/a33*a32*(2+d)/d*Dxx_p2;

    u_dx = Dx(u,dx); u_dx2 = Dx(u2,dx); u_dx1 = Dx(u1,dx);
    RHS = eps^2/(dt^2*a33^2)*p3 - 1/(dt*a33)*(d+2)/d*u_dx ...
        - 1/(dt*a33^2)*a31*(2+d)/d*u_dx1 - 1/(dt*a33^2)*a32*(2+d)/d*u_dx2 ...
        + a31/a33*(2+d)/d*Dxx(p1,dx) + a32/a33*(2+d)/d*Dxx(p2,dx);
    % Solve
    rhs = RHS'; 
    p_ = Lp3\rhs; p3 = p_'; 

  
    %% Compute u^n+1
    % p_dx3 = Lax_Friedrichs1D(p3,u,alpha_e,dx);
    % u3 = u  - dt*(a32*p_dx2 + a33*p_dx3)/eps;
    p_dx = Lax_Friedrichs1D(p1*a31 + p2*a32 + p3*a33,u,alpha_e,dx); % MIT LF!!!
    u3 = u - dt*p_dx;

    %% Compute rho^n+1
    u_dx_rho3 = Lax_Friedrichs1D(u3,rho,alpha_e,dx);
    rho3 = rho - dt/eps*(a32*u_dx_rho2 + a33*u_dx_rho3);
    % 
    % %% Compute T^n+1
    u_dx_T3 = Lax_Friedrichs1D(u3,T,alpha_e,dx);
    T3 = T3 - dt/eps*2/d*(a32*u_dx_T2 + a33*u_dx_T3);

    % rho3 = rho + d/(d+2)*(p3 - p3_star)*eps;
    % T3 = T3 + 2/d*d/(d+2)*(p3 - p3_star)*eps;

    %% UPDATE VARS
    T = T3; rho = rho3; u = u3; p = p3; 
    Gamma = Gamma3;
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
toc;