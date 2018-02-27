%% Homework 5, Part 2
clc; close all; clear all;
%% Create half annulus
nr = 10; nt = 30;
[r,theta] = meshgrid(linspace(0.5,1,nr),linspace(0,pi,nt));
x = r.*cos(theta); y = r.*sin(theta);
[X,Y] = meshgrid(x,y);

plot(x,y,'b-o',x',y','b-o'), axis equal

%% Iterate over elements on grid
% Note: We take the small elements on the rs grid as the new element on xy
%             Call the discretization of one element a subelement
A_tot = 0;
n_sub_elems = 10;
for i = 1:nr-1
    for j = 1:nt
        % Define corners
        ul_x = X(i,j);
        ur_x = X(i,j+1);
        ll_x = X(i+1,j);
        lr_x = X(i+1,j+1);
        xp = [ll_x lr_x ur_x ul_x];
        
        ul_y = Y(i,j);
        ur_y = Y(i,j+1);
        ll_y = Y(i+1,j);
        lr_y = Y(i+1,j+1);
        yp = [ll_y lr_y ur_y ul_y];

        % Pull discretization info
        [x_pts, y_pts, jac_arr, dA_arr] = rs2xy(xp,yp,n_sub_elems);
        
        % Add to total area
        A_tot = A_tot + sum(sum(dA_arr,1),2);
    end
end

disp(A_tot)