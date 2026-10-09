%% Computes divergence w. 4th order central diff.
% input: components of vector u = (u1,u2) \in R^2
% output: div u
function out = Div_O4(u1,u2,dx,dy)

    u1_h1 = zeros(size(u1,1),size(u1,2)+4);
    u1_h1(:,3:end-2) = u1; u1_h1(:,2) = u1(:,end-1); u1_h1(:,1) = u1(:,end-2);
    u1_h1(:,end-1) = u1(:,2); u1_h1(:,end) = u1(:,3);

    u2_h2 = zeros(size(u2,1)+4,size(u2,2));
    u2_h2(3:end-2,:) = u2; u2_h2(2,:) = u2(end-1,:); u2_h2(1,:) = u2(end-2,:);
    u2_h2(end-1,:) = u2(2,:); u2_h2(end,:) = u2(3,:);
    
    % Derivative in x-direction
    u1_dx = 1/(12*dx)*...
        (u1_h1(:,1:end-4)-8*u1_h1(:,2:end-3)+8*u1_h1(:,4:end-1)-u1_h1(:,5:end));
    % Derivative in y-direction
    u2_dy = 1/(12*dy)*...
        (u2_h2(1:end-4,:)-8*u2_h2(2:end-3,:)+8*u2_h2(4:end-1,:)-u2_h2(5:end,:));

    % Divergence
    out = u1_dx+u2_dy;

end