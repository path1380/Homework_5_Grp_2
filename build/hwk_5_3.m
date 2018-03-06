%% Homework 5, Part 2
clc; close all; clear all;

%% Study error
sample_grd_spacing = 50;
nt_arr = linspace(10,100,sample_grd_spacing);
error_arr = zeros(sample_grd_spacing);
for grd_spacing = 1:sample_grd_spacing

    %% Create half annulus
    nr = 10; %nt = 30;
    nt = nt_arr(grd_spacing);
    [r,theta] = meshgrid(linspace(0.5,1,nr),linspace(0,pi,nt));
    x = r.*cos(theta); y = r.*sin(theta);
    % [X,Y] = meshgrid(x,y);

    plot(x,y,'b-o',x',y','b-o'), axis equal
%     hold on
    %% Iterate over elements on grid
    % Note: We take the small elements on the rs grid as the new element on xy
    %             Call the discretization of one element a subelement
    A_tot = 0;
    n_sub_elems = 10;
    for i = 1:nt-1
        for j = 1:nr-1
            % Define corners
            ul_x = x(i,j);
            ur_x = x(i,j+1);
            ll_x = x(i+1,j);
            lr_x = x(i+1,j+1);
            xp = [ll_x lr_x ur_x ul_x];

            ul_y = y(i,j);
            ur_y = y(i,j+1);
            ll_y = y(i+1,j);
            lr_y = y(i+1,j+1);
            yp = [ll_y lr_y ur_y ul_y];

    %         % Visualize element being plotted
    %         scatter(xp, yp, 50, 'g', 'filled')

            % Pull discretization info
            [x_pts, y_pts, jac_arr, dA_arr] = rs2xy(xp,yp,n_sub_elems);

            % Add to total area
            A_tot = A_tot + sum(sum(dA_arr,1),2);
        end
    end

%     disp(A_tot)
    error_arr(grd_spacing) = 3/8*pi - A_tot;
end

figure
plot(nt_arr, error_arr)
hold on
plot(linspace(10,100), linspace(10,100).^(-2), '--')