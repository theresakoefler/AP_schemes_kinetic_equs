%% RK 3, WENO4
function [u1,u2,theta,X,Y,vort,dx,arr] = ShearLayer_RK3(eps,tau,n,tFinal)
    tic;
    delta = 0.05;
    rho = pi/15;
    CFL = 0.9;
    % Initialize Data Array
    arr = zeros(2000,2);

    % IMEX Parameter
    c = 1-1/sqrt(2); gamma = (c-1/2)/(c-1);

    a21_ = 1/2; 
    a31_ = 11/18; a32_ = 1/18;
    a41_ = 5/6; a42_ = -5/6; a43_ = 1/2;
    a51_ = 1/4; a52_ = 7/4; a53_ = 3/4; a54_ = -7/4;
        
    a22 = 1/2;
    a32 = 1/6; a33 = 1/2;
    a42 = -1/2; a43 = 1/2; a44 = 1/2;
    a52 = 3/2; a53 = -3/2; a54 = 1/2; a55 = 1/2;

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
    A2_(1,N-1) = 16; A2_(N,2) = 16; A2_(2,N-1) = -1; A2_(N,3) = -1; 
    A2_(1,N-2) = -1; A2_(N-1,2) =-1;

    % fourth order
    Lth4 = (kron(speye(N),sparse(A2_)) + kron(sparse(A2_),speye(N)))/(12*dx^2);
    % second order
    Lth2 = (kron(speye(N),sparse(A_)) + kron(sparse(A_),speye(N)))/dx^2;

    % Matrix for solving Helmholtz equ
    % internal step
    Lth_2 = (Lth4 ...
        - speye(N^2)*2*eps^2/(dt^2*a22^2));
    % internal step
    Lth_3 = (Lth4 - speye(N^2)*2*eps^2/(dt^2*a33^2));
    % internal step
    Lth_4 = (Lth4 - speye(N^2)*2*eps^2/(dt^2*a44^2));
    % final step
    Lth_5 = (Lth4 - speye(N^2)*2*eps^2/(dt^2*a55^2));
    
    while t<tFinal
        % correction of time step
        if t+dt>tFinal; dt=tFinal-t; end; t=t+dt; it = it+1;
        
        alpha_e = 1/100; 
        alpha_i = 0;

       %% ---------------- Internal Steps ---------------------------------
       %% ---------------- Step 1 -----------------------------------------
       v1_1 = v1; v2_1 = v2; u1_1 = u1; u2_1 = u2; theta_1 = theta; q_1 = q;

       %% ---------------- Step 2 -----------------------------------------
       %% compute v^(2)
       v1_flux_e = Lax_Friedrichs((-u1+4*eps^2*q),u2,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q),u2*0,v1,alpha_e,dx,dy);
       v1_2 = eps^2*tau/(eps^2*tau+dt*a22)*v1...
            - dt*tau*a21_/(4*(eps^2*tau+dt*c))*(v1_flux_e)...
            + dt*a21_/(eps^2*tau+dt*a22)*(1/2*(u2_1.^2-u1_1.^2));

       v2_flux_e = Lax_Friedrichs(u2,(u1+4*eps^2*q),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2*0,(q),v2,alpha_e,dx,dy);
       v2_2 = eps^2*tau/(eps^2*tau+dt*a22)*v2...
            - dt*tau*a21_/(4*(eps^2*tau+dt*a22))*(v2_flux_e)...
            + dt*a21_/(eps^2*tau+dt*a22)*(u1_1.*u2_1);

       %% Solve Helmholtz equation to compute theta^(2)
       div = Div_O4(u1,u2,dx,dy);
       % div = Lax_Friedrichs(u1,u2,theta,0,dx,dy);
       RHS = 2*eps^2/(dt^2*a22^2)*theta - 1/(dt*a22)*div...
        + Div_O4(Div_O4(-v1_2,v2_2,dx,dy),Div_O4(v2_2,v1_2,dx,dy),dx,dy);%DivDivAv_O4(v1_2,v2_2,dx,dy);

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
       [th_dx_2,th_dy_2] = Grad_O4(theta_2,dx,dy);
       flux_u1 = Lax_Friedrichs((-v1_2),v2_2,u1,alpha_e,dx,dy);
       u1_2 = u1 - dt*a22*(flux_u1) - dt*a22*th_dx_2;
    
       flux_u2 = Lax_Friedrichs(v2_2,(v1_2),u2,alpha_e,dx,dy);
       u2_2 = u2 - dt*a22*(flux_u2) - dt*a22*th_dy_2;

       %% compute q^n+1/2
       flux = Lax_Friedrichs(v1_2,v2_2,q,alpha_i,dx,dy);

       q_2 = eps^2*tau/(eps^2*tau+dt*a22)*q...
            - tau*dt*a22/(2*(eps^2*tau+dt*a22))*flux;


       %% --------------- Final Step ----------------------------------------
       %% compute v^n+1
       u1_ = (a31_*u1_1 + a32_* u1_2);
       u2_ = a31_*u2_1 + a32_* u2_2;
       q_ = a31_*q_1 + a32_* q_2;
       Fu1 = a31_*(1/2*(u2_1.^2-u1_1.^2))+a32_*(1/2*(u2_2.^2-u1_2.^2));
       Fu2 = a31_*u1_1.*u2_1+a32_*u1_2.*u2_2;

       v1_flux_e = Lax_Friedrichs((-u1_+4*eps^2*q_),u2_,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q_),u2_*0,v1,alpha_e,dx,dy);
       v1_3 = eps^2*tau/(eps^2*tau+dt*a33)*v1...
        - dt*tau/(4*(eps^2*tau+dt*a33))*(v1_flux_e)...
        - dt/(eps^2*tau+dt*a33)*a32*v1_2 ...
        + dt/(eps^2*tau+dt*a33)*(Fu1);

       v2_flux_e = Lax_Friedrichs(u2_,(u1_+4*eps^2*q_),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2_*0,(q_),v2,alpha_e,dx,dy);
       v2_3 = eps^2*tau/(eps^2*tau+dt*a33)*v2...
        - dt*tau/(4*(eps^2*tau+dt*a33))*(v2_flux_e)...
        - dt/(eps^2*tau+dt*a33)*a32*v2_2 ...
        + dt/(eps^2*tau+dt*a33)*Fu2;

       %% compute theta^n+1
       div_2 = Div_O4(u1_2,u2_2,dx,dy);

       % div = Div_O4(u1,u2,dx,dy);


       RHS = 2*eps^2/(dt^2*a33^2)*theta...
           - a32/(dt*a33^2)*div_2...
           - 1/(dt*a33)*div...
           + a33/a33*Div_O4(Div_O4(-v1_3,v2_3,dx,dy),Div_O4(v2_3,v1_3,dx,dy),dx,dy)...%DivDivAv_O4(v1_3,v2_3,dx,dy)...
           +a32/a33*Div_O4(Div_O4(-v1_2,v2_2,dx,dy),Div_O4(v2_2,v1_2,dx,dy),dx,dy)...%DivDivAv_O4(v1_2,v2_2,dx,dy)...
           +a32/a33*Lapl_O4(theta_2,dx,dy);




       rhs = reshape(-RHS,[],1);
       theta_ = Lth_3\rhs; theta_3 = reshape(theta_,N,N); 
    
       %% compute u^n+1
       v1_ = a32*v1_2 + a33*v1_3;
       v2_ = a32*v2_2 + a33*v2_3;
       % th_ = a32*theta_2 + a33*theta_3;

       flux_u1 = Lax_Friedrichs(-v1_,v2_,u1,alpha_e,dx,dy);
       flux_u2 = Lax_Friedrichs(v2_,v1_,u2,alpha_e,dx,dy);
       % [th_dx_2,th_dy_2] = Grad_O4(theta_2,dx,dy);
       [th_dx_3,th_dy_3] = Grad_O4(theta_3,dx,dy);

       u1_3 = u1 - dt*flux_u1-dt*(a32*th_dx_2 + a33*th_dx_3);
       u2_3 = u2 - dt*flux_u2-dt*(a32*th_dy_2 + a33*th_dy_3);

       %% compute q^n+1
       v1_ = a32*v1_2+a33*v1_3;
       v2_ = a32*v2_2+a33*v2_3;
       flux = Lax_Friedrichs(v1_,v2_,q,alpha_i,dx,dy);

       q_3 = eps^2*tau/(eps^2*tau+dt*a33)*q...
        - tau*dt/(2*(eps^2*tau+dt*a33))*(flux)...
        - dt/(eps^2*tau+dt*a33)*a32*q_2;

       %% FOURTH STEP
       %% compute v^n+1
       u1_ = (a41_*u1_1 + a42_* u1_2 + a43_* u1_3);
       u2_ = a41_*u2_1 + a42_* u2_2 + a43_* u2_3;
       q_ = a41_*q_1 + a42_* q_2 + a43_* q_3;
       Fu1 = a41_*(1/2*(u2_1.^2-u1_1.^2))+a42_*(1/2*(u2_2.^2-u1_2.^2))+a43_*(1/2*(u2_3.^2-u1_3.^2));
       Fu2 = a41_*u1_1.*u2_1+a42_*u1_2.*u2_2+a43_*u1_3.*u2_3;

       v1_flux_e = Lax_Friedrichs((-u1_+4*eps^2*q_),u2_,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q_),u2_*0,v1,alpha_e,dx,dy);
       v1_4 = eps^2*tau/(eps^2*tau+dt*a44)*v1...
        - dt*tau/(4*(eps^2*tau+dt*a44))*(v1_flux_e)...
        - dt/(eps^2*tau+dt*a44)*(a42*v1_2 + a43*v1_3) ...
        + dt/(eps^2*tau+dt*a44)*(Fu1);

       v2_flux_e = Lax_Friedrichs(u2_,(u1_+4*eps^2*q_),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2_*0,(q_),v2,alpha_e,dx,dy);
       v2_4 = eps^2*tau/(eps^2*tau+dt*a44)*v2...
        - dt*tau/(4*(eps^2*tau+dt*a44))*(v2_flux_e)...
        - dt/(eps^2*tau+dt*a44)*(a42*v2_2 + a43*v2_3) ...
        + dt/(eps^2*tau+dt*a44)*Fu2;

       %% compute theta^n+1
       % div_2 = Div_O4(u1_2,u2_2,dx,dy); 
       div_3 = Div_O4(u1_3,u2_3,dx,dy);
       % div_2 = Lax_Friedrichs(u1_2,u2_2,theta,0,dx,dy);

       % div = Div_O4(u1,u2,dx,dy);
   

       RHS = 2*eps^2/(dt^2*a44^2)*theta...
           - a42/(dt*a44^2)*div_2 - a43/(dt*a44^2)*div_3...
           - 1/(dt*a44)*div...
           + Div_O4(Div_O4(-v1_4,v2_4,dx,dy),Div_O4(v2_4,v1_4,dx,dy),dx,dy)...%DivDivAv_O4(v1_4,v2_4,dx,dy)...
           +a42/a44*Div_O4(Div_O4(-v1_2,v2_2,dx,dy),Div_O4(v2_2,v1_2,dx,dy),dx,dy)...%DivDivAv_O4(v1_2,v2_2,dx,dy)
           +a43/a44*Div_O4(Div_O4(-v1_3,v2_3,dx,dy),Div_O4(v2_3,v1_3,dx,dy),dx,dy)...%*DivDivAv_O4(v1_3,v2_3,dx,dy)...
           +a42/a44*Lapl_O4(theta_2,dx,dy)+a43/a44*Lapl_O4(theta_3,dx,dy);

       rhs = reshape(-RHS,[],1);
       theta_ = Lth_4\rhs; theta_4 = reshape(theta_,N,N); 
    
       %% compute u^n+1
       v1_ = a42*v1_2 + a43*v1_3 + a44*v1_4;
       v2_ = a42*v2_2 + a43*v2_3 + a44*v2_4;
       th_ = a42*theta_2 + a43*theta_3 + a44*theta_4;

       flux_u1 = Lax_Friedrichs(-v1_,v2_,u1,alpha_e,dx,dy);
       flux_u2 = Lax_Friedrichs(v2_,v1_,u2,alpha_e,dx,dy);
       % [th_dx_2,th_dy_2] = Grad_O4(theta_2,dx,dy);
       % [th_dx_3,th_dy_3] = Grad_O4(theta_3,dx,dy);
       [th_dx_4,th_dy_4] = Grad_O4(theta_4,dx,dy);

       u1_4 = u1 - dt*flux_u1-dt*(a42*th_dx_2 + a43*th_dx_3 + a44*th_dx_4);
       u2_4 = u2 - dt*flux_u2-dt*(a42*th_dy_2 + a43*th_dy_3 + a44*th_dy_4);

       %% compute q^n+1
       v1_ = a42*v1_2+a43*v1_3 +a44*v1_4;
       v2_ = a42*v2_2+a43*v2_3+a44*v2_4;
       flux = Lax_Friedrichs(v1_,v2_,q,alpha_i,dx,dy);

       q_4 = eps^2*tau/(eps^2*tau+dt*a44)*q...
        - tau*dt/(2*(eps^2*tau+dt*a44))*(flux)...
        - dt/(eps^2*tau+dt*a44)*(a42*q_2+a43*q_3);

       %% FIFTH STEP
       %% compute v^n+1
       u1_ = (a51_*u1_1 + a52_* u1_2 + a53_* u1_3 + a54_* u1_4);
       u2_ = a51_*u2_1 + a52_* u2_2 + a53_* u2_3 + a54_* u2_4;
       q_ = a51_*q_1 + a52_* q_2 + a53_* q_3 + a54_* q_4;
       Fu1 = a51_*(1/2*(u2_1.^2-u1_1.^2))+a52_*(1/2*(u2_2.^2-u1_2.^2))...
       +a53_*(1/2*(u2_3.^2-u1_3.^2))+a54_*(1/2*(u2_4.^2-u1_4.^2));
       Fu2 = a51_*u1_1.*u2_1+a52_*u1_2.*u2_2+a53_*u1_3.*u2_3+a54_*u1_4.*u2_4;

       v1_flux_e = Lax_Friedrichs((-u1_+4*eps^2*q_),u2_,v1,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs((q_),u2_*0,v1,alpha_e,dx,dy);
       v1_5 = eps^2*tau/(eps^2*tau+dt*a55)*v1...
        - dt*tau/(4*(eps^2*tau+dt*a55))*(v1_flux_e)...
        - dt/(eps^2*tau+dt*a55)*(a52*v1_2 + a53*v1_3 + a54*v1_4) ...
        + dt/(eps^2*tau+dt*a55)*(Fu1);

       v2_flux_e = Lax_Friedrichs(u2_,(u1_+4*eps^2*q_),v2,alpha_i,dx,dy)...
           +4*eps^2*Lax_Friedrichs(u2_*0,(q_),v2,alpha_e,dx,dy);
       v2_5 = eps^2*tau/(eps^2*tau+dt*a55)*v2...
        - dt*tau/(4*(eps^2*tau+dt*a55))*(v2_flux_e)...
        - dt/(eps^2*tau+dt*a55)*(a52*v2_2 + a53*v2_3 + a54*v2_4) ...
        + dt/(eps^2*tau+dt*a55)*Fu2;

       %% compute theta^n+1
       % div_2 = Div_O4(u1_2,u2_2,dx,dy); div_3 = Div_O4(u1_3,u2_3,dx,dy);
       div_4 = Div_O4(u1_4,u2_4,dx,dy);
       % div_2 = Lax_Friedrichs(u1_2,u2_2,theta,0,dx,dy);

       % div = Div_O4(u1,u2,dx,dy);
     

       RHS = 2*eps^2/(dt^2*a55^2)*theta...
           - a52/(dt*a55^2)*div_2 - a53/(dt*a55^2)*div_3 - a54/(dt*a55^2)*div_4...
           - 1/(dt*a55)*div...
           + Div_O4(Div_O4(-v1_5,v2_5,dx,dy),Div_O4(v2_5,v1_5,dx,dy),dx,dy)...%DivDivAv_O4(v1_5,v2_5,dx,dy)...
           +a52/a55*Div_O4(Div_O4(-v1_2,v2_2,dx,dy),Div_O4(v2_2,v1_2,dx,dy),dx,dy)...*DivDivAv_O4(v1_2,v2_2,dx,dy)
           +a53/a55*Div_O4(Div_O4(-v1_3,v2_3,dx,dy),Div_O4(v2_3,v1_3,dx,dy),dx,dy)...*DivDivAv_O4(v1_3,v2_3,dx,dy)...
           +a54/a55*Div_O4(Div_O4(-v1_4,v2_4,dx,dy),Div_O4(v2_4,v1_4,dx,dy),dx,dy)...*DivDivAv_O4(v1_4,v2_4,dx,dy)...
           +a52/a55*Lapl_O4(theta_2,dx,dy)+a53/a55*Lapl_O4(theta_3,dx,dy)...
           +a54/a55*Lapl_O4(theta_4,dx,dy);

       rhs = reshape(-RHS,[],1);
       theta_ = Lth_5\rhs; theta_5 = reshape(theta_,N,N); 
    
       %% compute u^n+1
       v1_ = a52*v1_2 + a53*v1_3 + a54*v1_4 + a55*v1_5;
       v2_ = a52*v2_2 + a53*v2_3 + a54*v2_4 + a55*v2_5;
       % th_ = a42*theta_2 + a43*theta_3 + a44*theta_4;

       flux_u1 = Lax_Friedrichs(-v1_,v2_,u1,alpha_e,dx,dy);
       flux_u2 = Lax_Friedrichs(v2_,v1_,u2,alpha_e,dx,dy);
       % [th_dx_2,th_dy_2] = Grad_O4(theta_2,dx,dy);
       % [th_dx_3,th_dy_3] = Grad_O4(theta_3,dx,dy);
       % [th_dx_4,th_dy_4] = Grad_O4(theta_4,dx,dy);
       [th_dx_5,th_dy_5] = Grad_O4(theta_5,dx,dy);

       u1_5 = u1 - dt*flux_u1-dt*(a52*th_dx_2 + a53*th_dx_3 + a54*th_dx_4 + a55*th_dx_5);
       u2_5 = u2 - dt*flux_u2-dt*(a52*th_dy_2 + a53*th_dy_3 + a54*th_dy_4 + a55*th_dy_5);

       %% compute q^n+1
       v1_ = a52*v1_2+a53*v1_3 +a54*v1_4 +a55*v1_5;
       v2_ = a52*v2_2+a53*v2_3+a54*v2_4 +a55*v2_5;
       flux = Lax_Friedrichs(v1_,v2_,q,alpha_i,dx,dy);

       q_5 = eps^2*tau/(eps^2*tau+dt*a55)*q...
        - tau*dt/(2*(eps^2*tau+dt*a55))*(flux)...
        - dt/(eps^2*tau+dt*a55)*(a52*q_2+a53*q_3+a54*q_4);

       %% UPDATE VARS
       u1 = u1_5; u2 = u2_5; v1 = v1_5; v2 = v2_5; theta = theta_5; q = q_5;
       %% Compute Divergence_______________________________________________
       % d = Div_O2(u1,u2,dx,dy);

       % arr(it+1,1) = t; arr(it+1,2) =  max(max(d));
        
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
        % contour(X,Y,vort,35)
        % title(['t is ',num2str(t)])
        % drawnow
    end
   
    %% SAVE DATA
    savename = ['rec/RK3_WENO5h_O4_tau=',int2str(tau),'eps=1e-6','_',int2str(n),];
    % savename = ['times/RK2WENO5_exp',int2str(anzahl)];
    save(savename, 'X','Y','u1','u2','dx','vort');

    % contour(X,Y,d,30)
    %     title(['t is ',num2str(t)])

toc;
end