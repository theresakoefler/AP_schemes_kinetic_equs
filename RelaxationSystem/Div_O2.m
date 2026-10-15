%% Computes divergence w. 2nd order central diff.
% input: components of vector u \in R^2
% output: div u
function out = Div_O2(u1,u2,dx,dy)
    u1_h1 = zeros(size(u1,1),size(u1,2)+2);
    u1_h1(:,2:end-1) = u1; u1_h1(:,1) = u1(:,end-1); u1_h1(:,end) = u1(:,2);

    u2_h2 = zeros(size(u2,1)+2,size(u2,2));
    u2_h2(2:end-1,:) = u2; u2_h2(1,:) = u2(end-1,:); u2_h2(end,:) = u2(2,:);
    
    u1_dx = 1/(2*dx)*(u1_h1(:,3:end)-u1_h1(:,1:end-2));
    u2_dy = 1/(2*dy)*(u2_h2(3:end,:)-u2_h2(1:end-2,:));

    out = u1_dx+u2_dy;

end