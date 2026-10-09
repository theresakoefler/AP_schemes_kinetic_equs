%% Compute Gradient w. 2nd order centr. diff
% input: theta \in R
% output: Gradient of theta [theta_x,theta_y]

function [th_dx,th_dy] = Grad_O2(theta,dx,dy)
    th_h1 = zeros(size(theta,1),size(theta,2)+2);
    th_h1(:,2:end-1) = theta; th_h1(:,1) = theta(:,end-1); 
    th_h1(:,end) = theta(:,2);

    th_h2 = zeros(size(theta,1)+2,size(theta,2));
    th_h2(2:end-1,:) = theta; th_h2(1,:) = theta(end-1,:); 
    th_h2(end,:) = theta(2,:);
    
    % Compute output
    th_dx = 1/(2*dx)*(th_h1(:,3:end)-th_h1(:,1:end-2));
    th_dy = 1/(2*dy)*(th_h2(3:end,:)-th_h2(1:end-2,:));

end