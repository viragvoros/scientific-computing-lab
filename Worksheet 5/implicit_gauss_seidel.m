% Function to solve 2D heat equation implicitly using Gauss-Seidel iterative solver

function [x,y,T] = implicit_gauss_seidel(L, n, tolerance, t_start, t_end, dt, current_T)

    x = linspace(0,L,n);        % x nodes
    y = linspace(0,L,n);        % y nodes
    dx = L/(n-1);               % gird size along x
    dy = L/(n-1);               % grid size along y
    n_t = (t_end-t_start) / dt;                 % number of time steps
    

    %Initialization
    if ~exist('current_T', 'var')
        T = ones(n, n);             % initializing T matrix
        T(:,1) = 0;                 % left boundary condition
        T(:,n) = 0;                 % right boundary condition
        T(1,:) = 0;                 % bottom boundary condition
        T(n,:) = 0;                 % top boundary condition
    else
        T = current_T;
    end
    
    T_old = T;                  % for updation old values in convergence loop
    T_prev = T;
    n_iteration = 1;            % to count the number of total iterations
    max_iteration = 1e5;        % maximum number of iterations

    k1 = dt/dx^2;               % for ease of calculation
    k2 = dt/dy^2;
    
    % time loop    
    for k = 1:n_t
        error = 1;              % error initialized to 1 before each time loop
        % convergence loop
        while error > tolerance && n_iteration < max_iteration
            % nodal loop
            for j = 2: (n-1)
                for i = 2: (n-1)
                    % implicit scheme using Gauss-Seidel
                    T(i,j) = (T_prev(i,j)+k1*(T(i-1,j)+T_old(i+1,j))+...
                        k2*(T(i,j-1)+T_old(i,j+1)))/(1+2*k1+2*k2);
                end
            end
            % convergence criterion
            error = max(max(abs(T - T_old)));
            % updating old values
            T_old = T;
            n_iteration = n_iteration +1;
        end
        % updating previous values
        T_prev = T;
    end
    
    if isinf(error) || isnan(error) || error > tolerance
        err.message = 'Could not find a solution';
        err.identifier = 'gauss_seidel_method:divergence';
        error(err);
    end
end