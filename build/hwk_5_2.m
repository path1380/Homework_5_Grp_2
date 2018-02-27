%% Homework 5, Part 2
clc; close all; clear all;
%% Create perturbed grid
n_x = 5;
n_y = n_x;

% Uniform grid
x = linspace(-1,1,n_x);
y = fliplr(linspace(-1,1,n_y));
[X,Y] = meshgrid(x,y);

% % Perturb grid
% Rp = R + rand(n_r,n_s)/25; 
% Sp = S + rand(n_r,n_s)/25;
% 
% % Clean up boundaries
% Rp(:,1) = -1;
% Rp(:,n_r) = 1;
% Sp(1,:) = 1;
% Sp(n_s,:) = -1;

% Temporarily keep grid uniform
Xp = X;
Yp = Y;

plot(Xp,Yp, 'b-x', Xp', Yp', 'b-x')

% Iterate over elements on grid
% Note: We take the small elements on the rs grid as the new element on xy
%             Call the discretization of one element a subelement
n_sub_elems = 5;
for i = 1:n_x-1
    for j = 1:n_y-1
        % Define corners
        ul_x = Xp(i,j);
        ur_x = Xp(i,j+1);
        ll_x = Xp(i+1,j);
        lr_x = Xp(i+1,j+1);
        xp = [ll_x lr_x ur_x ul_x];
        
        ul_y = Yp(i,j);
        ur_y = Yp(i,j+1);
        ll_y = Yp(i+1,j);
        lr_y = Yp(i+1,j+1);
        yp = [ll_y lr_y ur_y ul_y];

        % Pull discretization info
        [x_pts, y_pts, jac_arr, dA_arr] = rs2xy(xp,yp,n_sub_elems,2);
        
%         % Calculate Jacobian on each subelement
%         for m = 1:n_r-1
%             for n = 1:n_s-1
%                 delx = x_pts(1,m+1) - x_pts(1,m);
%                 dely = y_pts(m+1,n) - y_pts(m,n);
%                 delr = r_pts(1,m) - r_pts(1,m);
%                 dels = s_pts(m+1,n) - s_pts(m,n);
%                 
%                 J = [delx/delr delx/dels; dely/delr dely/dels];
%                 normJ = norm(J);
%                 
%             end
%         end
    end
end