%% Compute Laplacian with. 2nd order cen. diff
% input: theta \in R
% output: theta_xx + theta_yy

function out = Lapl_O2(theta,dx,dy)
    th_h1 = zeros(size(theta,1),size(theta,2)+2);
    th_h1(:,2:end-1) = theta; th_h1(:,1) = theta(:,end-1);
    th_h1(:,end) = theta(:,2);

    th_h2 = zeros(size(theta,1)+2,size(theta,2));
    th_h2(2:end-1,:) = theta; th_h2(1,:) = theta(end-1,:); 
    th_h2(end,:) = theta(2,:);

    th_xx = 1/dx^2*(th_h1(:,3:end)-2*theta+th_h1(:,1:end-2));
    th_yy = 1/dy^2*(th_h2(3:end,:)-2*theta+th_h2(1:end-2,:));

    out = th_xx + th_yy;
end