%% Homework 5, Part 2
clc; close all; clear all;
%% Create perturbed grid
n_r = 5;
n_s = n_r;

% Uniform grid
r = linspace(-1,1,n_r);
s = fliplr(linspace(-1,1,n_s));
[R,S] = meshgrid(r,s);

% Perturb grid
Rp = R + rand(n_r,n_s)/25; 
Sp = S + rand(n_r,n_s)/25;

% Clean up boundaries
Rp(:,1) = -1;
Rp(:,n_r) = 1;
Sp(1,:) = 1;
Sp(n_s,:) = -1;

plot(Rp,Sp, 'b-x', Rp', Sp', 'b-x')

% Iterate over elements on grid
% Note: We take the small elements on the rs grid as the new element on xy
%             Call the discretization of one element a subelement
n_sub_elems = 5;
for i = 1:n_r-1
    for j = 1:n_s-1
        % Define corners
        ul_x = Rp(i,j);
        ur_x = Rp(i,j+1);
        ll_x = Rp(i+1,j);
        lr_x = Rp(i+1,j+1);
        xp = [ll_x lr_x ur_x ul_x];
        
        ul_y = Sp(i,j);
        ur_y = Sp(i,j+1);
        ll_y = Sp(i+1,j);
        lr_y = Sp(i+1,j+1);
        yp = [ll_y lr_y ur_y ul_y];

        % Pull info from rs domain
        [x_pts, y_pts] = rs2xy(xp,yp,n_sub_elems,1);
        
    end
end