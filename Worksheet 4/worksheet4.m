close all
clear all
clc

N_size = [7 15 31 63 127];

% Running times
times = zeros(3, length(N_size));

% TODO implement storages
% Storages
storages = zeros(3, length(N_size));

% Error of Gauss-Seidel method
error = zeros(1, length(N_size));

for k = 1:length(N_size)
    
    Nx = N_size(k);
    Ny = N_size(k);
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);

    b = zeros(Nx, Ny);
    b_func = @(x,y) -2*pi^2 * sin(pi * x) * sin(pi * y);
    for i = 1:Nx
        for j = 1:Ny
        b(i,j) = b_func(i*hx, j*hy);
        end
    end

    % x(i,j) = sin(i*hx*pi) * sin(j*hy*pi)
    analytical = b ./  (-2*pi^2);

    x = zeros(Nx+2, Ny+2);

    % Full matrix time
    A = construct_matrix(Nx, Ny);
    b_flat = reshape(b, Nx*Ny, 1);
    tic;
    x(2:end-1, 2:end-1) = reshape(A\b_flat, Nx, Ny);
    times(1,k) = toc;

    % Sparse matrix
    A = sparse(A);
    tic;
    x(2:end-1, 2:end-1) = reshape(A\b_flat, Nx, Ny);
    times(2,k) = toc;
    
    % Gauss-Seidel time
    try
    tic;
    x(2:end-1, 2:end-1) = gauss_seidel_solve(Nx, Ny, b);
    times(3,k) = toc;
    catch err
        fprintf('Gauss-Seidel divergence at Nx %d Ny %d.', Nx, Ny);
        x = NaN;
    end

    % Calculting errors
    error(k) = sqrt(1/(Nx*Ny)*sum((x(2:end-1, 2:end-1)-analytical).^2, 'all'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the factor of error reduction if the step size is halved.
% Error reduction is the error of the previous and bigger time step divided
% by the given time step, hence an element-wise division of the array with
% itself but "shifted". A padding row is added because the first time step
% does not have any previous time step to compare the error with.

% The element at error_reduction[i] gives the error reduction from the
% i-1 to the i step.

error_reduction = [ 0, 	error(1:end-1) ./ error(2:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO
% Setting up the figures
% plot (contour or contourf ? COntourf seems prettier but WS says contour)
% contourf(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Displaying tables

full_matrix_table = table( ...
    [N_size(1:end-1); times(1, 1:end-1); storages(1, 1:end-1)], ...
    'VariableNames', {'direct solution with full matrix'}, ...
    'RowNames', {'Nx, Ny', 'runtime', 'storage'} ...
)

sparse_matrix_table = table( ...
    [N_size(1:end-1); times(2, 1:end-1); storages(2, 1:end-1)], ...
    'VariableNames', {'direct solution with sparse matrix'}, ...
    'RowNames', {'Nx, Ny', 'runtime', 'storage'} ...
)

gauss_seidel_table = table( ...
    [N_size(1:end-1); times(3, 1:end-1); storages(3, 1:end-1)], ...
    'VariableNames', {'iterative solution with Gauss-Seidel'}, ...
    'RowNames', {'Nx, Ny', 'runtime', 'storage'} ...
)

gauss_seidel_error_table = table( ...
    [N_size; error; error_reduction], ...
    'VariableNames', {'Errors of Gauss-Seidel method'}, ...
    'RowNames', {'Nx, Ny', 'error', 'error red.'} ...
)
