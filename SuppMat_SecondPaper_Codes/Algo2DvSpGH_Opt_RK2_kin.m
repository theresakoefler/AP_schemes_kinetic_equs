%% Algorithm 2D mit p - Trapezoidal Rule for Comp. Moments
% 2D-Test case, Mach Knudsen implicit, Diff exp, periodic BC
% ® Theresa Köfler
% Jan/2026
% Sec order in time, second order cen diff, third order weno

clear; close all; clc;
% function Algo2DvSpGH_Opt_RK2_kin(n,a)
eps = 1e-06;

tau = 0;

tFinal = 6;
CFL = 0.75;
vmax = 8;

d = 2; % dimension
a = 1; % LF-Parameter

n = 128; % number of grid pnts
M = 20; % # grid pnts velocity domain

% IMEX Verfahren - Typ A - erste Zeile d. Tableaus 0
c = 1-1/sqrt(2); gamma = (c-1/2)/(c-1);

a21_ = c; 
a31_ = 1-1/(2*c); a32_ = 1/(2*c);

a11 = 0;
a21 = 0; a22 = c;
a31 = 0; a32 = 1-gamma; a33 = gamma;

% Discretize spatial domain
x_min = 0; x_max = 2*pi; x = linspace(x_min,x_max,n); y = x; dx = x(2)-x(1); dy = dx;

% meshgrid
[X,Y] = meshgrid(x,y); 

% Discretize velocity space according to Gauss-Welsch-Algo.
beta = sqrt((1:M-1) / 2);
J = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(J);
[v1d, idx] = sort(diag(D)); % v_nodes
w1d = (V(1, idx).^2) * sqrt(pi);

% Erzeuge 2D-Gitter für Knoten
[V1, V2] = meshgrid(v1d, v1d);

% Erzeuge 2D-Gewichte (W_ij = w_i * w_j)
[W1, W2] = meshgrid(w1d, w1d);
W2D = W1 .* W2;
v1 = V1(:)*sqrt(2); % inkl skalierung f. Gauss-Hermit-Int
v2 = V2(:)*sqrt(2); % inkl skalierung f. Gauss-Hermit-Int
w  = W2D(:);


%% initial conditions
% initial condition like in Lattice-Boltzmann paper
delta = 0.05;
rho = pi/15;
u1 = (Y <= pi).*tanh(1/rho*(Y-pi/2)) + (Y> pi).*tanh(1/rho*(3*pi/2-Y));
u2 = delta*sin(X);

% initialize g
% for k=1:(M*M);
%      data.(strcat('f',num2str(k))) = zeros(n,n);
% end
M_sq = M*M;
F = zeros(n,n,M_sq); % als 3D array initialisieren
F1 = F; F2 = F; F3 = F;

% initialize p, rho, T
 rho = ones(n,n); T = rho*0; p = rho + T; S = T;


t = 0; dt = dx^2*CFL;%/vmax; 

% Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); A_(1,N-1) = 1; A_(N,2) = 1;
    r2 = zeros(1,N); r2(1) = -30; r2(2) = 16; r2(3) = -1; A2_ = toeplitz(r2); 
    A2_(1,N) = 16; A2_(N,1) = 16; A2_(2,N) = -1; A2_(N,2) = -1; A2_(1,N-1) = -1; A2_(N-1,1) =-1;
    LpH = (kron(speye(N),sparse(A2_)) + kron(sparse(A2_),speye(N)))/(12*dx^2);
    Lp = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

    % A_ = A_/dx^2;

% --- Before the while loop ---
mean_p_initial = mean(p(:));
while t < tFinal
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;

    alpha_e = a; alpha_i = 0;

    Lp2 = (speye(N^2)*eps^2/(dt^2*a22^2) - (d+2)/d*Lp);
    Lp3 = (speye(N^2)*eps^2/(dt^2*a33^2) - (d+2)/d*Lp);

    %% STAGE 1 ----------------------------------------------------
    rho_1 = rho; T_1 = T; p_1 = p; u1_1 = u1; u2_1 = u2; F1 = F;
    %% STAGE 2 ----------------------------------------------------
     % Wir formen F zu (n*n) x M_sq um für schnelle Matrix-Multiplikation
    F_reshaped1 = reshape(F1, n*n, M_sq); 
    % Gauss-Hermite int.
    coeff_base = (1/(2*pi)^(d/2)) * 2 * w;
    w_P2_11 = coeff_base .* (v1.*v1 - (v1.^2 + v2.^2)/2);
    w_P2_12 = coeff_base .* (v1.*v2);
    w_P2_21 = coeff_base .* (v1.*v2);
    w_P2_22 = coeff_base .* (v2.*v2 - (v1.^2 + v2.^2)/2);

    w_P3_1 = coeff_base.*v1.*(v1.^2+v2.^2 - (d+2));
    w_P3_2 = coeff_base.*v2.*(v1.^2+v2.^2 - (d+2));

    % Momente berechnen (Ergebnis ist n*n, dann zurück zu n x n)
    P1_11 = reshape(F_reshaped1 * w_P2_11, n, n);
    P1_12 = reshape(F_reshaped1 * w_P2_12, n, n);
    P1_21 = reshape(F_reshaped1 * w_P2_21, n, n);
    P1_22 = reshape(F_reshaped1 * w_P2_22, n, n);
    P1_1 = reshape(F_reshaped1 * w_P3_1, n, n);
    P1_2 = reshape(F_reshaped1 * w_P3_2, n, n);



     % terme berechnen
    rho_F1 = 1+eps*rho_1;
    u1_F1 = eps*u1_1./rho_F1; % Momentum!
    u2_F1 = eps*u2_1./rho_F1;
    T_F1 = 1 + eps*T_1./(1+eps*rho_1) - eps^2*(u1_1.^2+u2_1.^2)./(d*(1+eps*rho_1).^2);

    for k=1:(M*M);

    Avg_v1_1 = v1(k)*P1_11+v2(k)*P1_12;
    Avg_v2_1 = v1(k)*P1_21+v2(k)*P1_22;

    Bvg1_1 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P1_1;
    Bvg2_1 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P1_2;

    % M_101 = 1/(2*pi)^(d/2)*exp(-(v1(k)^2+v2(k)^2)/2);
    % M_F = rho_F./(2*pi*T_F).^(d/2).*exp(-((v1(k)-u1_F).^2+(v2(k)-u2_F).^2)./(2*T_F));

    Pf1 = rho_1 + v1(k)*u1_1+v2(k)*u2_1 + ((v1(k)^2+v2(k)^2)-d)/2*T_1;


    % Q = 1/eps^2*(M_F./M_101 - 1 - eps*Pf);
    %% -------
% 1. Vorfaktor berechnen (Dichte und Temperatur-Verhältnis)
% Beachte: (2*pi)^(d/2) kürzt sich fast weg, übrig bleibt:
    prefactor1 = rho_F1 ./ (T_F1.^(d/2));

% 2. Differenz der Exponenten (Vermeidung von Division durch M_101)
% Exp_F = -((v1(k)-u1_F).^2 + (v2(k)-u2_F).^2) ./ (2 * T_F)
% Exp_101 = -(v1(k)^2 + v2(k)^2) / 2
% Wir rechnen: Exp_F - Exp_101
    diff_exp1 = -((v1(k)-u1_F1).^2 + (v2(k)-u2_F1).^2) ./ (2 * T_F1) ...
           + (v1(k)^2 + v2(k)^2) / 2;

% 3. Das Verhältnis M_F / M_101 zusammenführen
    MF_over_M101_1 = prefactor1 .* exp(diff_exp1);

% 4. Den Q-Operator berechnen 
    if eps < 1e-04
    % Vorberechnungen für mehr Übersichtlichkeit
    v_sq = v1(k)^2 + v2(k)^2;         % |v|^2 (Skalar)
    u1_sq = u1_1.^2 + u2_1.^2;        % |u|^2 (n x n Matrix)
    vu1 = v1(k)*u1_1 + v2(k)*u2_1;    % v * u (n x n Matrix)
    
    % Taylor-Term für Q1 (entspricht exakt Q1 für eps -> 0)
    Q1 = 0.5 * vu1.^2 ...
          - 0.25 * v_sq * u1_sq ...
          + T_1 .* vu1 .* ((v_sq - 4) / 2) ...
          + (T_1.^2) .* (1 - v_sq + (v_sq^2)/8);

    else
    Q1 = 1/eps^2 * (MF_over_M101_1 - 1 - eps*Pf1);
    end
    %% -------


    Au1_1 = (v1(k)^2-(v1(k)^2+v2(k)^2)/d)*u1_1 + v1(k)*v2(k)*u2_1;
    Au2_1 = v1(k)*v2(k)*u1_1 + (v2(k)^2-(v1(k)^2+v2(k)^2)/d)*u2_1;
    BT1_1 = 1/2*v1(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_1;
    BT2_1 = 1/2*v2(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_1;


    % flux = Div_O2(Au1+BT1,Au2+BT2,dx,dy);

    flux1 = Div_O2(Au1_1+BT1_1+eps*(v1(k)*F1(:,:,k)- (Avg_v1_1+Bvg1_1)),...
        Au2_1+BT2_1 + eps*(v2(k)*F1(:,:,k) - (Avg_v2_1+Bvg2_1)),dx,dy);
    % flux = Lax_Friedrichs(Au1+BT1+eps*(v1(k)*F(:,:,k)- (Avg_v1+Bvg1)),...
    %     Au2+BT2 + eps*(v2(k)*F(:,:,k) - (Avg_v2+Bvg2)),F(:,:,k),1,dx,dy);

     F2(:,:,k) = ...
         eps^2*tau/(eps^2*tau+a22*dt)*F(:,:,k) ...
         - dt*tau/(eps^2*tau+a22*dt)*a21_*flux1 ...%- eps*dt*tau/(eps^2*tau+dt)*flux_f ...
         + dt/(eps^2*tau+a22*dt)*a21_*Q1;

    end

     % Wir formen F zu (n*n) x M_sq um für schnelle Matrix-Multiplikation
    F_reshaped2 = reshape(F2, n*n, M_sq); 

    % Momente berechnen (Ergebnis ist n*n, dann zurück zu n x n)

    P2_11 = reshape(F_reshaped2 * w_P2_11, n, n);
    P2_12 = reshape(F_reshaped2 * w_P2_12, n, n);
    P2_21 = reshape(F_reshaped2 * w_P2_21, n, n);
    P2_22 = reshape(F_reshaped2 * w_P2_22, n, n);
    P2_1 = reshape(F_reshaped2 * w_P3_1, n, n);
    P2_2 = reshape(F_reshaped2 * w_P3_2, n, n);
    
 

    %% Compute p^n+1/2
    flux_Gamma_2 = Div_O2(P2_1,P2_2,dx,dy);
    p_2 = p - dt/d*a22*flux_Gamma_2/eps;
    p_2star = p_2; % introduce helper var.
    T_2 = T - dt/d*a22*flux_Gamma_2;

    
    %% Compute u^n+1/2
    flux_u1_2 = Div_O2(P2_11,P2_12,dx,dy); % bei eps = 1, (tau = 0)
    flux_u2_2 = Div_O2(P2_21,P2_22,dx,dy);
    % 
    u1_2 = u1 -  dt*a22*flux_u1_2;
    u2_2 = u2 -  dt*a22*flux_u2_2;

    %% Compute p^n+1    
    RHS = eps^2/(dt^2*a22^2)*p_2 - 1/(dt*a22)*(d+2)/d*Div_O2(u1_2,u2_2,dx,dy);

    % % Solve
    % rhs = reshape(RHS,[],1); 
    % rho_ = Lp2\rhs; p_2 = reshape(rho_,N,N); 
    %% Stage 2 Solve
    rhs = reshape(RHS,[],1); 
    p_vec = Lp2\rhs; 
    p_2 = reshape(p_vec, N, N);% p_2_unpinned = p_2;
    
    % Differenz, BEVOR gepinnt wird
    % p_2star ist hier der "Predictor" vom Anfang der Stage
    % delta_p = p_2_unpinned - p_2star;
    % --- PIN THE MEAN HERE ---
    % p_2 = p_2 - mean(p_2(:)) + mean_p_initial;


    %% Compute u^n+1
    flux_u1_2LF = Lax_Friedrichs3(a22*P2_11,a22*P2_12,u1,alpha_e,dx,dy); 
    flux_u2_2LF = Lax_Friedrichs3(a22*P2_21,a22*P2_22,u2,alpha_e,dx,dy);
    [p_dx2,p_dy2] = Grad_O2(p_2,dx,dy);

    u1_2 = u1 -  dt*flux_u1_2LF  - dt*a22*p_dx2;
    u2_2 = u2 -  dt*flux_u2_2LF  - dt*a22*p_dy2;

    %% Compute T^n+1 and rho^n+1 algebraically
    rho_2 = rho + d/(d+2)*(eps*p_2 - eps*p_2star);
    T_2 = T_2 + 2/d*d/(d+2)*(eps*p_2 - eps*p_2star);

    % %% Compute T^n+1
    % T_2 = T_2 - dt*2/d*Lax_Friedrichs1(a22*u1_2/eps,a22*u2_2/eps,T,alpha_e*0.65,dx,dy);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute rho^n+1
    % rho_2 = rho - dt*Lax_Friedrichs1(a22*u1_2/eps,a22*u2_2/eps,rho,alpha_e*0.65,dx,dy);%Div_O2(u1,u2,dx,dy);


    %% STAGE 3 ----------------------------------------------------

     % terme berechnen
    rho_F2 = 1+eps*rho_2;
    u1_F2 = eps*u1_2./rho_F2; % Momentum!
    u2_F2 = eps*u2_2./rho_F2;
    T_F2 = 1 + eps*T_2./(1+eps*rho_2) - eps^2*(u1_2.^2+u2_2.^2)./(d*(1+eps*rho_2).^2);

    for k=1:(M*M);

    Avg_v1_1 = v1(k)*P1_11+v2(k)*P1_12;
    Avg_v2_1 = v1(k)*P1_21+v2(k)*P1_22;
    Avg_v1_2 = v1(k)*P2_11+v2(k)*P2_12;
    Avg_v2_2 = v1(k)*P2_21+v2(k)*P2_22;

    Bvg1_1 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P1_1;
    Bvg2_1 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P1_2;
    Bvg1_2 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P2_1;
    Bvg2_2 = ((v1(k)^2+v2(k)^2)-d)/(2*d)*P2_2;


    Pf1 = rho_1 + v1(k)*u1_1+v2(k)*u2_1 + ((v1(k)^2+v2(k)^2)-d)/2*T_1;
    Pf2 = rho_2 + v1(k)*u1_2+v2(k)*u2_2 + ((v1(k)^2+v2(k)^2)-d)/2*T_2;


    % Q = 1/eps^2*(M_F./M_101 - 1 - eps*Pf);
    %% -------
% 1. Vorfaktor berechnen (Dichte und Temperatur-Verhältnis)
% Beachte: (2*pi)^(d/2) kürzt sich fast weg, übrig bleibt:
    prefactor1 = rho_F1 ./ (T_F1.^(d/2));
    prefactor2 = rho_F2 ./ (T_F2.^(d/2));

% 2. Differenz der Exponenten (Vermeidung von Division durch M_101)
% Exp_F = -((v1(k)-u1_F).^2 + (v2(k)-u2_F).^2) ./ (2 * T_F)
% Exp_101 = -(v1(k)^2 + v2(k)^2) / 2
% Wir rechnen: Exp_F - Exp_101
    diff_exp1 = -((v1(k)-u1_F1).^2 + (v2(k)-u2_F1).^2) ./ (2 * T_F1) ...
           + (v1(k)^2 + v2(k)^2) / 2;
    diff_exp2 = -((v1(k)-u1_F2).^2 + (v2(k)-u2_F2).^2) ./ (2 * T_F2) ...
           + (v1(k)^2 + v2(k)^2) / 2;

% 3. Das Verhältnis M_F / M_101 zusammenführen
    MF_over_M101_1 = prefactor1 .* exp(diff_exp1);
    MF_over_M101_2 = prefactor2 .* exp(diff_exp2);

% 4. Den Q-Operator berechnen
    if eps < 1e-04
    % Vorberechnungen für mehr Übersichtlichkeit
    v_sq = v1(k)^2 + v2(k)^2;         % |v|^2 (Skalar)
    u1_sq = u1_1.^2 + u2_1.^2;        % |u|^2 (n x n Matrix)
    vu1 = v1(k)*u1_1 + v2(k)*u2_1;    % v * u (n x n Matrix)
    
    % Taylor-Term für Q1 (entspricht exakt Q1 für eps -> 0)
    Q1 = 0.5 * vu1.^2 ...
              - 0.25 * v_sq * u1_sq ...
              + T_1 .* vu1 .* ((v_sq - 4) / 2) ...
              + (T_1.^2) .* (1 - v_sq + (v_sq^2)/8);
    
    % Analog für Stage 3 (Q2):
    u2_sq = u1_2.^2 + u2_2.^2;
    vu2 = v1(k)*u1_2 + v2(k)*u2_2;
    Q2 = 0.5 * vu2.^2 ...
              - 0.25 * v_sq * u2_sq ...
              + T_2 .* vu2 .* ((v_sq - 4) / 2) ...
              + (T_2.^2) .* (1 - v_sq + (v_sq^2)/8);
    else

    Q1 = 1/eps^2 * (MF_over_M101_1 - 1 - eps*Pf1);
    Q2 = 1/eps^2 * (MF_over_M101_2 - 1 - eps*Pf2);
    end
    %% -------


    Au1_1 = (v1(k)^2-(v1(k)^2+v2(k)^2)/d)*u1_1 + v1(k)*v2(k)*u2_1;
    Au2_1 = v1(k)*v2(k)*u1_1 + (v2(k)^2-(v1(k)^2+v2(k)^2)/d)*u2_1;
    BT1_1 = 1/2*v1(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_1;
    BT2_1 = 1/2*v2(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_1;

    Au1_2 = (v1(k)^2-(v1(k)^2+v2(k)^2)/d)*u1_2 + v1(k)*v2(k)*u2_2;
    Au2_2 = v1(k)*v2(k)*u1_2 + (v2(k)^2-(v1(k)^2+v2(k)^2)/d)*u2_2;
    BT1_2 = 1/2*v1(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_2;
    BT2_2 = 1/2*v2(k)*((v1(k)^2+v2(k)^2) - (d+2))*T_2;


    % flux = Div_O2(Au1+BT1,Au2+BT2,dx,dy);

    flux1 = Div_O2(Au1_1+BT1_1+eps*(v1(k)*F1(:,:,k)- (Avg_v1_1+Bvg1_1)),...
        Au2_1+BT2_1 + eps*(v2(k)*F1(:,:,k) - (Avg_v2_1+Bvg2_1)),dx,dy);
    flux2 = Div_O2(Au1_2+BT1_2+eps*(v1(k)*F2(:,:,k)- (Avg_v1_2+Bvg1_2)),...
        Au2_2+BT2_2 + eps*(v2(k)*F2(:,:,k) - (Avg_v2_2+Bvg2_2)),dx,dy);
    % flux = Lax_Friedrichs(Au1+BT1+eps*(v1(k)*F(:,:,k)- (Avg_v1+Bvg1)),...
    %     Au2+BT2 + eps*(v2(k)*F(:,:,k) - (Avg_v2+Bvg2)),F(:,:,k),1,dx,dy);

     F3(:,:,k) = ...
         eps^2*tau/(eps^2*tau+a33*dt)*F(:,:,k) ...
         - dt*tau/(eps^2*tau+a33*dt)*(a31_*flux1 + a32_*flux2) ...%- eps*dt*tau/(eps^2*tau+dt)*flux_f ...
         + dt/(eps^2*tau+a33*dt)*(a31_*Q1 + a32_*Q2) ...
         - dt/(eps^2*tau+a33*dt)*a32*F2(:,:,k);

    end

     % Wir formen F zu (n*n) x M_sq um für schnelle Matrix-Multiplikation
    F_reshaped3 = reshape(F3, n*n, M_sq); 

    % Momente berechnen (Ergebnis ist n*n, dann zurück zu n x n)

    P3_11 = reshape(F_reshaped3 * w_P2_11, n, n);
    P3_12 = reshape(F_reshaped3 * w_P2_12, n, n);
    P3_21 = reshape(F_reshaped3 * w_P2_21, n, n);
    P3_22 = reshape(F_reshaped3 * w_P2_22, n, n);
    P3_1 = reshape(F_reshaped3 * w_P3_1, n, n);
    P3_2 = reshape(F_reshaped3 * w_P3_2, n, n);



    %% Compute p^n+1/2
        flux_Gamma_3 = Div_O2(P3_1,P3_2,dx,dy);
    p_3 = p - dt/d*(a32*flux_Gamma_2 + a33*flux_Gamma_3)/eps;
    p_3star = p_3;
    T_3 = T - dt/d*(a32*flux_Gamma_2 + a33*flux_Gamma_3);


    %% Compute u^n+1/2
    flux_u1_3 = Div_O2(P3_11,P3_12,dx,dy); % bei eps = 1, (tau = 0)
    flux_u2_3 = Div_O2(P3_21,P3_22,dx,dy);
    % 
    u1_3 = u1 -  dt*(a32*flux_u1_2 + a33*flux_u1_3);
    u2_3 = u2 -  dt*(a32*flux_u2_2 + a33*flux_u2_3);

    %% Compute p^n+1    
    RHS = eps^2/(dt^2*a33^2)*p_3 - 1/(dt*a33)*(d+2)/d*Div_O2(u1_3,u2_3,dx,dy)...
        - 1/(dt*a33^2)*a32*(d+2)/d*Div_O2(u1_2,u2_2,dx,dy)...
        + a32/a33*(d+2)/d*Lapl_O2(p_2,dx,dy);

    % % Solve
    % rhs = reshape(RHS,[],1); 
    % rho_ = Lp3\rhs; p_3 = reshape(rho_,N,N); 
    %% Stage 3 Solve
    rhs = reshape(RHS,[],1); 
    p_vec = Lp3\rhs; 
    p_3 = reshape(p_vec, N, N); %p_3_unpinned = p_3;
    
    % Differenz, BEVOR gepinnt wird
    % p_2star ist hier der "Predictor" vom Anfang der Stage
    %delta_p = p_3_unpinned - p_3star;
    
    % --- PIN THE MEAN HERE ---
    %p_3 = p_3 - mean(p_3(:)) + mean_p_initial;


    %% Compute u^n+1
    comp11 = a32*P2_11 + a33*P3_11;
    comp12 = a32*P2_12 + a33*P3_12;
    comp21 = a32*P2_21 + a33*P3_21;
    comp22 = a32*P2_22 + a33*P3_22;
    flux_u1_3LF = Lax_Friedrichs3(comp11,comp12,u1,alpha_e,dx,dy);
    flux_u2_3LF = Lax_Friedrichs3(comp21,comp22,u2,alpha_e,dx,dy);
    [p_dx3,p_dy3] = Grad_O2(p_3,dx,dy);

    u1_3 = u1 - dt*flux_u1_3LF  - dt*(a32*p_dx2 + a33*p_dx3);
    u2_3 = u2 - dt*flux_u2_3LF   - dt*(a32*p_dy2 + a33*p_dy3);

    %% Compute rho^n+1 and T^n+1 by using algebraic properties
    rho_3 = rho + d/(d+2)*(eps*p_3 - eps*p_3star);
    T_3 = T_3 + 2/d*d/(d+2)*(eps*p_3 - eps*p_3star);

    % %% Compute T^n+1
    % T_3 = T_3 - dt*2/d*Lax_Friedrichs1((a32*u1_2+a33*u1_3)/eps,(a32*u2_2+a33*u2_3)/eps,T,alpha_e*0.65,dx,dy);%Div_O2(u1,u2,dx,dy);
    % 
    % %% Compute rho^n+1
    % rho_3 = rho - dt*Lax_Friedrichs1((a32*u1_2+a33*u1_3)/eps,(a32*u2_2+a33*u2_3)/eps,rho,alpha_e*0.65,dx,dy);%Div_O2(u1,u2,dx,dy);

    %% ________UPDATE VARS_______________________________
    rho = rho_3; T = T_3; p = p_3; u1 = u1_3; u2 = u2_3; F = F3;
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

    contour(X,Y,vort,34)
    title(['t is ',num2str(t)])
    drawnow
end
% %% SAVE DATA
%     folder = 'Data_2D_ConvTest_2nd_Paper';
%     if ~exist(folder, 'dir')
%         mkdir(folder);
%     end
% 
%     savename = ['Data_2D_ConvTest_2nd_Paper/RK2WENO3O2LF_kin',int2str(a),'_tau=0','eps=1e-6','_',int2str(n),'.mat'];
%     % savename = ['1DTests-IMEX_RK2/RK2kWENO3O2_tau=01','eps=1','M=',int2str(M),'n=',int2str(n)];
%     save(savename, 'X','Y','dx','rho','u1','u2','T','p');
% end