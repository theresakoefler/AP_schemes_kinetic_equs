%% Algorithm with v-Space in 1D, IMEX RK2, GH-Quadratur 
% 1D-Test case, Mach Knudsen implicit, Diff exp
% ® Theresa Köfler
% Nov/2025


clear; close all; clc;
% function Algo1DvSpGHRK2(n)
tic;
eps = 1e-0;

tau = 0.01;

tFinal = 0.1;
CFL = 0.75;
vmax = 10;

d = 1; % dimension

n = 500; % number of grid pnts
M = 20; % # grid pnts velocity domain

% IMEX Parameter
% c = 2.25; w2_ = 1/(2*c);
% 
% a21_ = c/w2_; 
% a31_ = 1-w2_; a32_ = w2_;
% 
% a11 = c;
% a21 = 0; a22 = c;
% a31 = 1-c; a32 = 0; a33 = c;

% Anderes IMEX Verfahren - Typ A - erste Zeile d. Tableaus 0
c = 1-1/sqrt(2); gamma = (c-1/2)/(c-1);

a21_ = c; 
a31_ = 1-1/(2*c); a32_ = 1/(2*c);

a11 = 0;
a21 = 0; a22 = c;
a31 = 0; a32 = 1-gamma; a33 = gamma;

% Discretize spatial domain
x_min = 0; x_max = 1; x = linspace(x_min,x_max,n); dx = x(2)-x(1);

% Discretize velocity space according to Gauss-Welsch-Algo.
beta = sqrt((1:M-1) / 2);
J = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(J);
[v1d, idx] = sort(diag(D)); % v_nodes
w1d = (V(1, idx).^2) * sqrt(pi);


v = v1d*sqrt(2); % inkl skalierung f. Gauss-Hermit-Int
w  = w1d(:);



% v = linspace(-vmax, vmax, M);
% dv = v(2)-v(1);
% w1d = ones(1, M) * dv;
% % w1d(1) = dv/2; w1d(end) = dv/2; % 1D Trapez-Korrektur


%% initial conditions
u = x*0;

% initialize F
F = zeros(M,n);
F1 = zeros(M,n); F2 = F1; F3 = F1;
% initialize rho, 
% rho = ones(1,n);
% 
% rho(0 <= x & x <= 0.2) = 1;
% rho(0.2 < x & x <= 0.3) = 1+eps^2;
% rho(0.3 < x & x <= 0.7) = 1;
% rho(0.7 < x & x <= 0.8) = 1-eps^2;
% % % rho = sin(2*pi*x);
% % % initialize T
% T = u*0; p = (rho+T)/eps; %Bg = T;

% Andere IC
rho = zeros(1,n); T = rho; p = rho+T;
for i=1:n
    rho(i)= sin(x(i)*2*pi)+0.5; 
end
p = (rho+T)/eps; u = T;

t = 0; dt = dx*CFL*dx; 

% Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); A_(1,N-1) = 1; A_(N,2) = 1;
    r2 = zeros(1,N); r2(1) = -30; r2(2) = 16; r2(3) = -1; A2_ = toeplitz(r2); 
    A2_(1,N) = 16; A2_(N,1) = 16; A2_(2,N) = -1; A2_(N,2) = -1; A2_(1,N-1) = -1; A2_(N-1,1) =-1;
    % LpH = (kron(speye(N),sparse(A2_)) + kron(sparse(A2_),speye(N)))/(12*dx^2);
    % Lp = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

    A_ = A_/dx^2;
while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    
    alpha_e = 1; alpha_i = 0;

 
    F1 = F; u1 = u; rho1 = rho; p1 = p; T1 = T;

    %% Compute g^n+1
    %% Gauss-Hermite
    Maxw = 1/(2*pi)^(d/2)*exp(-(abs(v.^2)/2));
    poly_v = v.*(abs(v).^2 - (d+2)); % M x 1 vector
    P3_1 = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F1; % 1 x n vector
    % P3_1 = (Maxw.*w.*poly_v)'*F1; % 1 x n vector -> Riemann!!

    V_mat = repmat(v, 1, size(u,2)); % [M x n]
    % reusable vars for Q-calc
    rho_F1 = 1+eps*rho1;
    u_F1 = eps*u1./rho_F1; % Momentum!
    T_F1 = 1 + eps*T1./(1+eps*rho1) - eps^2*abs(u1).^2./(d*(1+eps*rho1).^2);

    Bg = 0*u; % initialize

    Bvg1 = (abs(v).^2-d)/(2*d)*P3_1; % Part of ((I-P)(vg)) decomp.

    num_diff =  alpha_e;
    flux_f = Lax_Friedrichs1D(a21_*(V_mat.*F1 - Bvg1), ...
        F,num_diff,dx); % mit LF!
    % flux_f = Dx(a21_*(V_mat.*F1 - Bvg1),dx); % ohne LF!


    % LBM Coll Op
    % B = (abs(v(k))^2-(d+2))*T*v(k).*u;
    % Q = 1/(2)*( B); % LBM type Q 
    if eps < 1e-3
        % Term 3: T * v * u * (v^2 - 3)/2
        term3 = T1 .* u1 .* (V_mat .* (V_mat.^2 - 3) / 2);
        
        % Term 4: T^2 * (v^4/8 - 3/4 * v^2 + 3/8)
        term4 = (T1.^2) .* (V_mat.^4 / 8 - 0.75 * V_mat.^2 + 3/8);

        % Gesamtes Q1
        Q1_taylor = term3 + term4;
        Q1 = Q1_taylor;
    else
    % Anderes Q
        M_101 = 1/(2*pi)^(d/2)*exp(-(abs(V_mat).^2)/2);
        M_F1 = rho_F1./(2*pi*T_F1).^(d/2).*exp(-(abs(V_mat-u_F1).^2)./(2*T_F1));
    
        Pf1 = rho1 + V_mat.*u1 + (abs(V_mat).^2-d)/2.*T1;    
        Q1 = 1/eps^2*(M_F1./M_101 - 1 - eps*Pf1);
    end

    % flux calculation
    BT1 = 1/2*v.*(abs(v).^2 - (d+2))*T1; % v.*(abs(v).^2 - (d+2))*T/2 -> M x n Matrix

    flux = Dx(a21_*BT1,dx);

     F2 = ...
         eps^2*tau/(eps^2*tau+dt*a22)*F ...
         - dt*tau/(eps^2*tau+dt*a22)*flux - eps*dt*tau/(eps^2*tau+dt*a22)*flux_f ...
         + dt/(eps^2*tau+dt*a22)*a21_*Q1;


      % Maxw = 1/(2*pi)^(d/2)*exp(-(abs(v(k))^2)/2);
      Bg = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F2;
      % Bg = (Maxw.*w.*poly_v)'*F2;


    %% Compute p^n+1/2, T^n+1/2
    T2 =T - dt*a22*Dx(Bg,dx);
    p2 = p - dt/d*a22*Dx(Bg,dx)/eps;
    p2_star = p2;


    %% Compute p^n+1
    Lp2 = (speye(N)*eps^2/(dt^2*a22^2) - (d+2)/d*A_);
    u_dx = Dx(u,dx);
    RHS = eps^2/(dt^2*a22^2)*p2 - 1/(dt*a22)*(d+2)/d*u_dx;

    % Solve
    rhs = RHS'; 
    p_ = Lp2\rhs; p2 = p_'; 


    %% Compute u^n+1
    p_dx = Lax_Friedrichs1D(p2*a22,u,alpha_e,dx); % MIT LF!!!
    % p_dx = Dx(p,dx);

    u2 = u - dt*p_dx;


    rho2 = rho + d/(d+2)*(p2 - p2_star)*eps;
    T2 = T2 + 2/d*d/(d+2)*(p2 - p2_star)*eps;


    %% Compute g^n+1
    %% Gauss-Hermite

    poly_v = v.*(abs(v).^2 - (d+2)); % M x 1 vector
    P3_2 = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F2; % 1 x n vector
    % P3_2 = (Maxw.*w.*poly_v)'*F2; % 1 x n vector

    % reusable vars for Q-calc
    rho_F2 = 1+eps*rho2;
    u_F2 = eps*u2./rho_F2; % Momentum!
    T_F2 = 1 + eps*T2./(1+eps*rho2) - eps^2*abs(u2).^2./(d*(1+eps*rho2).^2);

    Bg2 = 0*u; Bg3 = Bg2; Bg1 = Bg3;% initialize
   

    Bvg1 = (abs(v).^2-d)/(2*d)*P3_1; % Part of ((I-P)(vg)) decomp.
    Bvg2 = (abs(v).^2-d)/(2*d)*P3_2; % Part of ((I-P)(vg)) decomp.

    num_diff =  alpha_e;
    flux_f = Lax_Friedrichs1D(a31_*(V_mat.*F1 - Bvg1)+a32_*(V_mat.*F2 - Bvg2), ...
        F,num_diff,dx); % mit LF!
    % flux_f = Dx(a31_*(V_mat.*F1 - Bvg1)+a32_*(V_mat.*F2 - Bvg2),dx); % ohne LF!

    % LBM Coll Op
    % B = (abs(v(k))^2-(d+2))*T*v(k).*u;
    % Q = 1/(2)*( B); % LBM type Q 
    
    if eps < 1e-3
        % Term 3: T * v * u * (v^2 - 3)/2
        term3_1 = T1 .* u1 .* (V_mat .* (V_mat.^2 - 3) / 2);
        term3_2 = T2 .* u2 .* (V_mat .* (V_mat.^2 - 3) / 2);
        
        % Term 4: T^2 * (v^4/8 - 3/4 * v^2 + 3/8)
        term4_1 = (T1.^2) .* (V_mat.^4 / 8 - 0.75 * V_mat.^2 + 3/8);
        term4_2 = (T2.^2) .* (V_mat.^4 / 8 - 0.75 * V_mat.^2 + 3/8);

        % Gesamtes Q1
        Q1_taylor = term3_1 + term4_1;
        Q2_taylor = term3_2 + term4_2;
        Q1 = Q1_taylor;
        Q2 = Q2_taylor;
    else
        % Anderes Q
        M_101 = 1/(2*pi)^(d/2)*exp(-(abs(V_mat).^2)/2);
        M_F1 = rho_F1./(2*pi*T_F1).^(d/2).*exp(-(abs(V_mat-u_F1).^2)./(2*T_F1));
        M_F2 = rho_F2./(2*pi*T_F2).^(d/2).*exp(-(abs(V_mat-u_F2).^2)./(2*T_F2));
    
        Pf1 = rho1 + V_mat.*u1 + (abs(V_mat).^2-d)/2.*T1;  
        Pf2 = rho2 + V_mat.*u2 + (abs(V_mat).^2-d)/2.*T2;
        Q1 = 1/eps^2*(M_F1./M_101 - 1 - eps*Pf1);
        Q2 = 1/eps^2*(M_F2./M_101 - 1 - eps*Pf2);
    end

    % flux calculation
    BT1 = 1/2*v.*(abs(v).^2 - (d+2))*T1; % v.*(abs(v).^2 - (d+2))*T/2 -> M x n Matrix
    BT2 = 1/2*v.*(abs(v).^2 - (d+2))*T2;
    flux = Dx(a31_*BT1 + a32_*BT2,dx);

     F3 = ...
         eps^2*tau/(eps^2*tau+dt*a33)*F ...
         - dt*tau/(eps^2*tau+dt*a33)*flux - eps*dt*tau/(eps^2*tau+dt*a33)*flux_f ...
         + dt/(eps^2*tau+dt*a33)*(a31_*Q1 + a32_*Q2)...
         - dt/(eps^2*tau+dt*a33)*(a32*F2+a31*F1);


      % Maxw = 1/(2*pi)^(d/2)*exp(-(abs(v(k))^2)/2);
      Bg1 = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F1;
      Bg2 = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F2;
      Bg3 = 1/(2*pi)^(d/2)*sqrt(2)*(w.*poly_v)'*F3;

      % Bg1 = (Maxw.*w.*poly_v)'*F1;
      % Bg2 = (Maxw.*w.*poly_v)'*F2;
      % Bg3 = (Maxw.*w.*poly_v)'*F3;



    %% Compute p^n+1/2, T^n+1/2

    T3 =T - dt*(a31*Dx(Bg1,dx)+ a32*Dx(Bg2,dx) + a33*Dx(Bg3,dx));
    p3 = p - dt/d*(a31*Dx(Bg1,dx)+a32*Dx(Bg2,dx) + a33*Dx(Bg3,dx))/eps;
    p3_star = p3;

    %% Compute p^n+1
    Lp3 = (speye(N)*eps^2/(dt^2*a33^2) - (d+2)/d*A_);
    u_dx = Dx(u,dx); u_dx2 = Dx(u2,dx); u_dx1 = Dx(u1,dx);
     RHS = eps^2/(dt^2*a33^2)*p3 - 1/(dt*a33)*(d+2)/d*u_dx ...
        - 1/(dt*a33^2)*a31*(2+d)/d*u_dx1 - 1/(dt*a33^2)*a32*(2+d)/d*u_dx2 ...
        + a31/a33*(2+d)/d*Dxx(p1,dx) + a32/a33*(2+d)/d*Dxx(p2,dx);

    % % Solve
    rhs = RHS'; 
    p_ = Lp3\rhs; p3 = p_'; 




    %% Compute u^n+1
    p_dx = Lax_Friedrichs1D(p1*a31 + p2*a32 + p3*a33,u,alpha_e,dx); % MIT LF!!!
    % p_dx = Dx(p,dx);

    u3 = u - dt*p_dx;


    rho3 = rho + d/(d+2)*(p3 - p3_star)*eps;
    T3 = T3 + 2/d*d/(d+2)*(p3 - p3_star)*eps;


    F = F3;
    rho = rho3; T = T3; 
    u = u3; p =  p3; 
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

    % subplot(4,1,4)
    % plot(x,Bg3);
    % title('Bg3')

    sgtitle(['t is ',num2str(t)])
    drawnow
end

% %% SAVE DATA
%     savename = ['Data_1D_ConvTest_2nd_Paper/RK2WENO3O2_tau=0','eps=001','_',int2str(n),'t01'];
%     % savename = ['1DTests-IMEX_RK2/RK2kWENO3O2_tau=01','eps=1','M=',int2str(M),'n=',int2str(n)];
%     save(savename, 'x','rho','u','T','p','dx');
%  toc;
% end