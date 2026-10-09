%% Compute Laplacian with. 4th order cen. diff
% input: theta \in R
% output: theta_xx + theta_yy

function out = Lapl_O4(theta,dx,dy)
    th_h1 = zeros(size(theta,1),size(theta,2)+4);
    th_h1(:,3:end-2) = theta; th_h1(:,2) = theta(:,end-1); 
    th_h1(:,1) = theta(:,end-2);
    th_h1(:,end-1) = theta(:,2); th_h1(:,end) = theta(:,3);

    th_h2 = zeros(size(theta,1)+4,size(theta,2));
    th_h2(3:end-2,:) = theta; th_h2(2,:) = theta(end-1,:); 
    th_h2(1,:) = theta(end-2,:);
    th_h2(end-1,:) = theta(2,:); th_h2(end,:) = theta(3,:);

    % Second order derivatives
    th_xx = 1/(12*dx^2)*...
        (-th_h1(:,1:end-4)+16*th_h1(:,2:end-3)-30*theta...
        +16*th_h1(:,4:end-1)-th_h1(:,5:end));
    th_yy = 1/(12*dy^2)*...
        (-th_h2(1:end-4,:)+16*th_h2(2:end-3,:)-30*theta...
        +16*th_h2(4:end-1,:)-th_h2(5:end,:));

    out = th_xx + th_yy;
end