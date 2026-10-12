% function [u1,u2,X,Y,vort,dx,arr] = ShearLayer_RK2(eps,tau,n,tFinal)
%     tic;
clear; eps = 1e-6; tau = 0; n = 180; tFinal = 6;
    delta = 0.05;
    rho = pi/50;
    CFL = 0.75;
    % Initialize Data Array
    arr = zeros(2000,2);

    % IMEX Parameter
    c = 2.25; gamma = (c-1/2)/(c-1);

    a21_ = c; 
    a31_ = 1-1/(2*c); a32_ = 1/(2*c);
        
    a22 = c;
    a32 = 1-gamma; a33 = gamma;

    % Discretize spatial domain
    x_min = 0; x_max = 2*pi; 
    x = linspace(x_min,x_max,n); y = x; dx = x(2)-x(1); dy = dx;

    % meshgrid
    [X,Y] = meshgrid(x,y); 

    % initial conditions
    u1 = (Y <= pi).*tanh(1/rho*(Y-pi/2)) + (Y> pi).*tanh(1/rho*(3*pi/2-Y));
    u2 = delta*sin(X);

    theta = u1*0; v1 = u1*0; v2 = u1*0; q = u1*0;

    %% ______ SOLVER LOOP _________________________________________________
    
    t = 0; it = 0; dt = dx^2*CFL;
    % 2D Poisson Matrix
    N = n; 
    r = zeros(1,N); r(1) = -2; r(2) = 1; A_ = toeplitz(r); 
    A_(1,N-1) = 1; A_(N,2) = 1;
    r2 = zeros(1,N); r2(1) = -30; r2(2) = 16; r2(3) = -1; A2_ = toeplitz(r2); 
    A2_(1,N) = 16; A2_(N,1) = 16; A2_(2,N) = -1; A2_(N,2) = -1; 
    A2_(1,N-1) = -1; A2_(N-1,1) =-1;

    % fourth order
    Lth4 = (kron(speye(N),sparse(A2_)) + kron(sparse(A2_),speye(N)))/(12*dx^2);
    % second order
    Lth2 = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

    % Matrix for solving Helmholtz equ
    % internal step
    Lth_2 = (Lth2 ...
        - speye(N^2)*2*eps^2/(dt^2*a22^2));
    % final step
    Lth = (Lth2 - speye(N^2)*2*eps^2/(dt^2*a33^2));
    
    while t<tFinal
        % correction of time step
        if t+dt>tFinal; dt=tFinal-t; end; t=t+dt; it = it+1;
        
        alpha_e = 1; 
        alpha_i = 0;

       %% ---------------- Internal Steps ---------------------------------
       %% ---------------- Step 1 -----------------------------------------
       v1_1 = v1; v2_1 = v2; u1_1 = u1; u2_1 = u2; theta_1 = theta; q_1 = q;

       %% ---------------- Step 2 -----------------------------------------
       %% compute v^(2)
       v1_flux_e = Lax_Friedrichs((-u1+4*eps^2*q),u2,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q),u2*0,v1,alpha_e,dx,dy);
       v1_2 = eps^2*tau/(eps^2*tau+dt*c)*v1...
            - dt*tau*c/(4*(eps^2*tau+dt*c))*(v1_flux_e)...
            + dt*c/(eps^2*tau+dt*c)*(1/2*(u2_1.^2-u1_1.^2));

       v2_flux_e = Lax_Friedrichs(u2,(u1+4*eps^2*q),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2*0,(q),v2,alpha_e,dx,dy);
       v2_2 = eps^2*tau/(eps^2*tau+dt*c)*v2...
            - dt*tau*c/(4*(eps^2*tau+dt*c))*(v2_flux_e)...
            + dt*c/(eps^2*tau+dt*c)*(u1_1.*u2_1);

       %% Solve Helmholtz equation to compute theta^(2)
       div = Div_O2(u1,u2,dx,dy);
       % div = Lax_Friedrichs(u1,u2,theta,0,dx,dy);
       RHS = 2*eps^2/(dt^2*c^2)*theta - 1/(dt*c)*div...
        +DivDivAv_O2(v1_2,v2_2,dx,dy);

       %______
       % flux_u1 = Lax_Friedrichs((-v1_2),v2_2,u1,alpha_e,dx,dy);
       % u1_ = u1 - dt*c*(flux_u1);% - dt*c*th_dx_2;
       % 
       % flux_u2 = Lax_Friedrichs(v2_2,(v1_2),u2,alpha_e,dx,dy);
       % u2_ = u2 - dt*c*(flux_u2);% - dt*c*th_dy_2;
       % 
       % flux = Lax_Friedrichs(u1_,u2_,theta,0,dx,dy);
       % RHS = 2*eps^2/(dt^2*c^2)*theta - 1/(dt*c)*flux;

       rhs = reshape(-RHS,[],1);
       theta_ = Lth_2\rhs; theta_2 = reshape(theta_,N,N); 

       %% compute u^(2)
       [th_dx_2,th_dy_2] = Grad_O2(theta_2,dx,dy);
       flux_u1 = Lax_Friedrichs((-v1_2),v2_2,u1,alpha_e,dx,dy);
       u1_2 = u1 - dt*c*(flux_u1) - dt*c*th_dx_2;
    
       flux_u2 = Lax_Friedrichs(v2_2,(v1_2),u2,alpha_e,dx,dy);
       u2_2 = u2 - dt*c*(flux_u2) - dt*c*th_dy_2;

       %% compute q^n+1/2
       flux = Lax_Friedrichs(v1_2,v2_2,q,alpha_i,dx,dy);

       q_2 = eps^2*tau/(eps^2*tau+dt*c)*q...
            - tau*dt*c/(2*(eps^2*tau+dt*c))*flux;


       %% --------------- Final Step ----------------------------------------
       %% compute v^n+1
       u1_ = (a31_*u1_1 + a32_* u1_2);
       u2_ = a31_*u2_1 + a32_* u2_2;
       q_ = a31_*q_1 + a32_* q_2;
       Fu1 = a31_*(1/2*(u2_1.^2-u1_1.^2))+a32_*(1/2*(u2_2.^2-u1_2.^2));
       Fu2 = a31_*u1_1.*u2_1+a32_*u1_2.*u2_2;

       v1_flux_e = Lax_Friedrichs((-u1_+4*eps^2*q_),u2_,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q_),u2_*0,v1,alpha_e,dx,dy);
       v1 = eps^2*tau/(eps^2*tau+dt*gamma)*v1...
        - dt*tau/(4*(eps^2*tau+dt*gamma))*(v1_flux_e)...
        - dt/(eps^2*tau+dt*gamma)*(1-gamma)*v1_2 ...
        + dt/(eps^2*tau+dt*gamma)*(Fu1);

       v2_flux_e = Lax_Friedrichs(u2_,(u1_+4*eps^2*q_),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2_*0,(q_),v2,alpha_e,dx,dy);
       v2 = eps^2*tau/(eps^2*tau+dt*gamma)*v2...
        - dt*tau/(4*(eps^2*tau+dt*gamma))*(v2_flux_e)...
        - dt/(eps^2*tau+dt*gamma)*(1-gamma)*v2_2 ...
        + dt/(eps^2*tau+dt*gamma)*Fu2;

       %% compute theta^n+1
       div_2 = Div_O2(u1_2,u2_2,dx,dy);
       % div_2 = Lax_Friedrichs(u1_2,u2_2,theta,0,dx,dy);

       div = Div_O2(u1,u2,dx,dy);
       % div = Lax_Friedrichs(u1,u2,theta,0,dx,dy);

       RHS = 2*eps^2/(dt^2*gamma^2)*theta...
           - (1-gamma)/(dt*gamma^2)*div_2...
           - 1/(dt*gamma)*div...
           + DivDivAv_O2(v1,v2,dx,dy)...
           +(1-gamma)/(gamma)*DivDivAv_O2(v1_2,v2_2,dx,dy)...
           +(1-gamma)/(gamma)*Lapl_O2(theta_2,dx,dy);

       %_________
       % v1_ = (1-gamma)*v1_2 + gamma*v1;
       % v2_ = (1-gamma)*v2_2 + gamma*v2;
       % %th_ = (1-gamma)*theta_2 + gamma*theta;
       % 
       % flux_u1 = Lax_Friedrichs(-v1_,v2_,u1,alpha_e,dx,dy);
       % flux_u2 = Lax_Friedrichs(v2_,v1_,u2,alpha_e,dx,dy);
       % 
       % u1_ = u1 - dt*flux_u1;%-dt*((1-gamma)*th_dx_2 + gamma*th_dx);
       % u2_ = u2 - dt*flux_u2;%-dt*((1-gamma)*th_dy_2 + gamma*th_dy);
       % 
       % flux = Lax_Friedrichs((1-gamma)*u1_2+gamma*u1_,...
       %     (1-gamma)*u2_2+gamma*u2_,theta,0,dx,dy);
       % RHS = 2*eps^2/(dt^2*gamma^2)*theta...
       %     - 1/(dt*gamma^2)*flux...
       %     +(1-gamma)/(gamma)*Lapl_O4(theta_2,dx,dy);


       rhs = reshape(-RHS,[],1);
       theta_ = Lth\rhs; theta = reshape(theta_,N,N); 
    
       %% compute u^n+1
       v1_ = (1-gamma)*v1_2 + gamma*v1;
       v2_ = (1-gamma)*v2_2 + gamma*v2;
       th_ = (1-gamma)*theta_2 + gamma*theta;

       flux_u1 = Lax_Friedrichs(-v1_,v2_,u1,alpha_e,dx,dy);
       flux_u2 = Lax_Friedrichs(v2_,v1_,u2,alpha_e,dx,dy);
       [th_dx_2,th_dy_2] = Grad_O2(theta_2,dx,dy);
       [th_dx,th_dy] = Grad_O2(theta,dx,dy);

       u1 = u1 - dt*flux_u1-dt*((1-gamma)*th_dx_2 + gamma*th_dx);
       u2 = u2 - dt*flux_u2-dt*((1-gamma)*th_dy_2 + gamma*th_dy);

       %% compute q^n+1
       v1_ = (1-gamma)*v1_2+gamma*v1;
       v2_ = (1-gamma)*v2_2+gamma*v2;
       flux = Lax_Friedrichs(v1_,v2_,q,alpha_i,dx,dy);

       q = eps^2*tau/(eps^2*tau+dt*gamma)*q...
        - tau*dt*gamma/(2*(eps^2*tau+dt*gamma))*(flux)...
        - dt/(eps^2*tau+dt*gamma)*(1-gamma)*q_2;

       %% Compute Divergence_______________________________________________
       d = Div_O2(u1,u2,dx,dy);

       arr(it+1,1) = t; arr(it+1,2) =  max(max(d));
        
       %% ___________ PLOT ________________________________________________
        % du1_dy = 1/(2*dx)*(u1(3:end,2:end-1) - u1(1:end-2,2:end-1));
        % du2_dx = 1/(2*dx)*(u2(2:end-1,3:end) - u2(2:end-1,1:end-2));
        % contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),du2_dx - du1_dy,35)
        % title(['t is ',num2str(t)])
        % 
        % drawnow
        %% compute Vort
        u2_(:,2:n+1) = u2;
        u2_(:,1) = u2(:,end-1);
        u2_(:,n+2) = u2(:,2);
        
        u1_(2:n+1,:) = u1;
        u1_(1,:) = u1(end-1,:);
        u1_(n+2,:) = u1(2,:);
        du1_dy = 1/(2*dx)*(u1_(3:end,:) - u1_(1:end-2,:));
        du2_dx = 1/(2*dx)*(u2_(:,3:end) - u2_(:,1:end-2));
        vort = du2_dx - du1_dy;

        %% Plot 
        contour(X,Y,vort,35)
        title(['t is ',num2str(t)])
        drawnow
    end
   
%     %% SAVE DATA
%     savename = ['rec3/RK2WENO3O2_tau=',int2str(tau),'eps=1e-6','_',int2str(n),];
%     % savename = ['times/RK2WENO5_exp',int2str(anzahl)];
%     save(savename, 'X','Y','u1','u2','dx','vort');
% 
%     % contour(X,Y,d,30)
%     %     title(['t is ',num2str(t)])
% 
% toc;
% end


