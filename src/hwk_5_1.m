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

% % corner 1 = 1, corners 2 = 3 = 4 = 0 (lower left)
% ll = [1 0 0 0]';
% kvec1 = inv(A)*ll;
% 
% % corner 2 = 1, corners 1 = 3 = 4 = 0 (lower right)
% lr = [0 1 0 0]';
% kvec2 = inv(A)*lr;
% 
% % corner 3
% ur = [0 0 1 0]';
% kvec3 = inv(A)*ur;
% 
% % corner 4
% ul = [0 0 0 1]';
% kvec4 = inv(A)*ul;

%% Map reference domain to x-y plane
% x and y coordinates
x = [-2 5 4 -1]';
y = [-4 -2 2 5]';

% r and s coordinates
r = [-1 1 1 -1];
s = [-1 -1 1 1];

% % subplot(1, 2, 1);
% figure
% plot(r(1:2), s(1:2), 'k'), axis equal
% hold on
% plot(r(2:3), s(2:3), 'bl'), axis equal
% plot(r(3:4), s(3:4), 'g'), axis equal
% plot((r(4) r(1)), (s(4) s(1)), 'r'), axis equal

% plot(r(2:3), s(2:3), 'bl'), axis equal


