% Function to solve 2D heat equation explicitly

function [x,y,T] = explicit(L, n, t_start, t_end, dt, current_T)

    x = linspace(0,L,n);        % x nodes
    y = linspace(0,L,n);        % y nodes
    dx = L/(n-1);               % gird size along x
    dy = L/(n-1);               % grid size along y
    n_t = (t_end-t_start) / dt;             % number of time steps


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

    k1 = dt/dx^2;               % for ease of calculation
    k2 = dt/dy^2;

    % time loop    
    for k = 1:n_t   
        % nodal loop
        for j = 2: (n-1)
            for i = 2: (n-1)
                % explicit scheme for 2D heat equation
                T(i,j) = (1-2*k1-2*k2) * T_old(i,j) + k1*(T_old(i+1,j) + T_old(i-1,j))...
                    + k2*(T_old(i,j+1) + T_old(i,j-1)); 
            end
        end
        % updating old values
        T_old = T;

    end
end