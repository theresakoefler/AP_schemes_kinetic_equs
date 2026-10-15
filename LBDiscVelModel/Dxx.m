function th_dx = Dxx(theta,dx)
    th_h1 = zeros(size(theta,1),size(theta,2)+2);
    th_h1(:,2:end-1) = theta; th_h1(:,1) = theta(:,end-1); 
    th_h1(:,end) = theta(:,2);

    
    % Compute output
    th_dx = 1/(dx^2)*(th_h1(:,3:end)-2*th_h1(:,2:end-1)+th_h1(:,1:end-2));

end