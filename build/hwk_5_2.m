%% Homework 5, Part 2
clc; close all; clear all;
%% Create perturbed grid
n_x = 10;
n_y = n_x;

% Uniform grid
x = linspace(-1,1,n_x);
y = fliplr(linspace(-1,1,n_y));
[X,Y] = meshgrid(x,y);

% Perturb grid
Xp = X + rand(n_x,n_y)/25; 
Yp = Y + rand(n_x,n_y)/25;

% Clean up boundaries
Xp(:,1) = -1;
Xp(:,n_x) = 1;
Yp(1,:) = 1;
Yp(n_y,:) = -1;

plot(Xp,Yp, 'b-x', Xp', Yp', 'b-x')

% Iterate over elements on grid
% Note: We take the small elements on the rs grid as the new element on xy
%             Call the discretization of one element a subelement
A_tot = 0;
n_sub_elems = 3;
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
        [x_pts, y_pts, jac_arr, dA_arr] = rs2xy(xp,yp,n_sub_elems);
        
        % Add to total area
        A_tot = A_tot + sum(sum(dA_arr,1),2);
    end
end

disp(A_tot)