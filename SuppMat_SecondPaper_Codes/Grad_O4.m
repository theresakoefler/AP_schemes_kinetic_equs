%% Compute Gradient w. 4th order centr. diff
% input: theta \in R
% output: Gradient of theta [theta_x,theta_y]

function [th_dx,th_dy] = Grad_O4(theta,dx,dy)
    th_h1 = zeros(size(theta,1),size(theta,2)+4);
    th_h1(:,3:end-2) = theta; th_h1(:,2) = theta(:,end-1); 
    th_h1(:,1) = theta(:,end-2);
    th_h1(:,end-1) = theta(:,2); th_h1(:,end) = theta(:,3);

    th_h2 = zeros(size(theta,1)+4,size(theta,2));
    th_h2(3:end-2,:) = theta; th_h2(2,:) = theta(end-1,:); 
    th_h2(1,:) = theta(end-2,:);
    th_h2(end-1,:) = theta(2,:); th_h2(end,:) = theta(3,:);
    
    % Derivative in x-direction
    th_dx = 1/(12*dx)*...
        (th_h1(:,1:end-4)-8*th_h1(:,2:end-3)+8*th_h1(:,4:end-1)-th_h1(:,5:end));
    % Derivative in y-direction
    th_dy = 1/(12*dy)*...
        (th_h2(1:end-4,:)-8*th_h2(2:end-3,:)+8*th_h2(4:end-1,:)-th_h2(5:end,:));

end