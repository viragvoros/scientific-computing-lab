% Solving 2D heat equation explicitly and implicitly with Gauss-Seidel

clear all
close all
clc

L = 1;                      % domain length along x and y                                          
n = 3;                      % number of grid points along x and y
tolerance = 1e-6;           % convergence criterion
t = 0.125;                  % desired time of temperature distribution
dt_explicit = 1/64;         % time step size for exlicit method
dt_imp_gauss = 1/64;        % time step size for implicit gauss-seidel

explicit(L, n, t, dt_explicit); %function for solving using explicit scheme
implicit_gauss_seidel(L, n, tolerance, t, dt_imp_gauss); % implicit gauss-seidel

