%% Homework 5, Part 1
clc; close all; clear all;

%% Bilinear mapping
% Find a bilinear mapping from (r,s) e [-1,1]^2 into a straight sided quad
% Assume bilinear eqn of the form
%   f(r,s) = k1 + k2*r + k3*s + k4*r*s
%   f(r,s) = A*k

% 4 _ _ _ 3
%  |     |
%  |     |
%  |_ _ _|
% 1       2

% Find A by evaluating f(r,s) at all corners
A = [1 -1 -1 1; 1 1 -1 -1; 1 1 1 1; 1 -1 1 -1];

%% Map reference domain to x-y plane
% x and y coordinates of the quadrilateral
xp = [-2 5 4 -1]';
yp = [-4 -2 2 5]';

x_coeffs = inv(A)*xp;
y_coeffs = inv(A)*yp;

% Discretize rs square
n_points = 5;
r12 = linspace(-1,1,n_points);
s12 = linspace(-1,-1,n_points);
r23 = linspace(1,1,n_points);
s23 = linspace(-1,1,n_points);
r34 = linspace(1,-1,n_points);
s34 = linspace(1,1,n_points);
r41 = linspace(-1,-1,n_points);
s41 = linspace(1,-1,n_points);

% Discretize xy quadrilateral
x12 = zeros(numel(r12), 1);
y12 = zeros(numel(r12), 1);
x23 = zeros(numel(r23), 1);
y23 = zeros(numel(r23), 1);
x34 = zeros(numel(r34), 1);
y34 = zeros(numel(r34), 1);
x41 = zeros(numel(r41), 1);
y41 = zeros(numel(r41), 1);

% Side 1-2
for i = 1:numel(r12)
    x12(i) = x_coeffs(1) + x_coeffs(2)*r12(i) + x_coeffs(3)*s12(i) ...
                        + x_coeffs(4)*r12(i)*s12(i);
    y12(i) = y_coeffs(1) + y_coeffs(2)*r12(i) + y_coeffs(3)*s12(i) ...
                        + y_coeffs(4)*r12(i)*s12(i);
end

% Side 2-3
for i = 1:numel(r23)
    x23(i) = x_coeffs(1) + x_coeffs(2)*r23(i) + x_coeffs(3)*s23(i) ...
                        + x_coeffs(4)*r23(i)*s23(i);
    y23(i) = y_coeffs(1) + y_coeffs(2)*r23(i) + y_coeffs(3)*s23(i) ...
                        + y_coeffs(4)*r23(i)*s23(i);
end

% Side 3-4
for i = 1:numel(r34)
    x34(i) = x_coeffs(1) + x_coeffs(2)*r34(i) + x_coeffs(3)*s34(i) ...
                        + x_coeffs(4)*r34(i)*s34(i);
    y34(i) = y_coeffs(1) + y_coeffs(2)*r34(i) + y_coeffs(3)*s34(i) ...
                        + y_coeffs(4)*r34(i)*s34(i);
end

% Side 4-1
for i = 1:numel(r41)
    x41(i) = x_coeffs(1) + x_coeffs(2)*r41(i) + x_coeffs(3)*s41(i) ...
                        + x_coeffs(4)*r41(i)*s41(i);
    y41(i) = y_coeffs(1) + y_coeffs(2)*r41(i) + y_coeffs(3)*s41(i) ...
                        + y_coeffs(4)*r41(i)*s41(i);
end

%% Plot
% Format plot
figure = gcf;
figure.PaperUnits = 'inches';
figure.PaperPosition = [0 0 12 3];

% r-s square
subplot(1,2,1)
plot(r12, s12, 'o'), axis equal
hold on
plot(r23, s23, 'o'), axis equal
plot(r34, s34, 'o'), axis equal
plot(r41, s41, 'o'), axis equal

% x-y quad
subplot(1,2,2)
plot(x12, y12, 'o'), axis equal
hold on
plot(x23, y23, 'o'), axis equal
plot(x34, y34, 'o'), axis equal
plot(x41, y41, 'o'), axis equal

saveas(gcf, 'foobar', 'png')