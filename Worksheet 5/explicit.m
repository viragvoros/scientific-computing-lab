% Function to solve 2D heat equation explicitly

function [x, y, T] = explicit(L, n, t_start, t_end, dt, current_T)

    x = linspace(0, L, n+2);      % x nodes
    y = linspace(0, L, n+2);      % y nodes
    dx = L / (n+1);               % grid size along x
    dy = L / (n+1);               % grid size along y
    n_t = (t_end - t_start) / dt; % number of time steps


    % Initialization
    if ~exist('current_T', 'var')
        T = ones(n+2, n+2);       % initializing T matrix
        T(:, 1) = 0;              % left boundary condition
        T(:, end) = 0;            % right boundary condition
        T(1, :) = 0;              % bottom boundary condition
        T(end, :) = 0;            % top boundary condition
    else
        T = current_T;
    end
        
    T_old = T;                    % updation old values in convergence loop

    k1 = dt / dx^2;               % for ease of calculation
    k2 = dt / dy^2;

    % time loop    
    for k = 1:n_t   
        T(2:end-1, 2:end-1) = (1 - 2*k1 - 2*k2) * T(2:end-1, 2:end-1) ...
            + k1 * (T(1:end-2, 2:end-1) + T(3:end, 2:end-1)) ...
            + k2 * (T(2:end-1, 1:end-2) + T(2:end-1, 3:end));
    end
end
