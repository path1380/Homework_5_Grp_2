%% Homework 5, Part 2
clc; close all; clear all;
%% Create perturbed grid
n_r = 10;
n_s = n_r;

% Uniform grid
r = linspace(-1,1,n_r);
s = linspace(-1,1,n_s);
[R,S] = meshgrid(r,s);

% Perturb grid
Rp = R + rand(n_r,n_s)/25; 
Sp = S + rand(n_r,n_s)/25;

% Clean up boundaries
Rp(:,1) = -1;
Rp(:,n_r) = 1;
Sp(1,:) = -1;
Sp(n_s,:) = 1;

plot(Rp,Sp, 'b-x', Rp', Sp', 'b-x')