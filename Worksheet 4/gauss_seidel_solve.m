function x = gauss_seidel_solve(Nx, Ny, b)

% Solve the system Ax = b on the unit square with homogenous Dirichlet BC. A is the matrix representing
% the laplacian operator with finite differences on a rectangular grid with Nx inner nodes on the horizontal axis
% and Ny inner nodes on the vertical axis.
% The Gauss-Seidel algorithm is used, but the matrix is never explicitly
% built.
% Inputs:
%   Nx, Ny : dimensions of the grid (number of inner nodes on both axes)
%   b : Nx * Ny matrix with the right hand side value at a given point of
%   the grid
% Outputs: 
%   x : Nx * Ny matrix with the nodal values in grid format

nb_dofs = Nx * Ny;
x = zeros(nb_dofs, 1); % The output will get reshaped to a Nx * Ny matrix in the end
b = reshape(b, nb_dofs, 1); % Convert the Nx x Ny matrix to a (Nx*Ny) column vector

hx = 1/(Nx+1);
hy = 1/(Ny+1);

% "Stencil" of the laplacian (coefficients in the matrix)
main_term = -2/hx^2 + -2/hy^2;
horizontal_term = 1/hx^2;
vertical_term = 1/hy^2;

% Standard gauss_seidel iteration is of the form :
% x_i^(n+1)  = 1/A(i,i) * (b(i) - sum(j < i) A(i,j)*x_j^(n+1) - sum(j > i)
% A(i,j)*x_j^n
% The "n+1" terms on the right hand side are automatically used if each DOF
% is updated in increasing order
residual = norm(b);
tol_residual = 1e-4;
max_iter = 1e5;
n = 0;

while residual > tol_residual && n < max_iter
    % Updating x
    for i = 1:Nx
        for j = 1:Ny
            rhs = 0;
            if i ~= 1
                rhs = rhs + horizontal_term * x(i-1 + (j-1) * Nx);
            end
            if i ~= Nx
                rhs = rhs + horizontal_term * x(i+1 + (j-1) * Nx);
            end
            if j ~= 1
                rhs = rhs + vertical_term * x(i + (j-2) * Nx);
            end
            if j ~= Ny
                rhs = rhs + vertical_term * x(i + (j) * Nx);
            end
            x(i + (j-1) * Nx) = 1/main_term * (b(i + (j-1) * Nx) - rhs);
        end
    end
    
    %Recomputing residual = |Ax - b|
    residual_vec = -b;
    
    for i = 1:Nx
        for j = 1:Ny
            rhs = 0;
            if i ~= 1
                rhs = rhs + horizontal_term * x(i-1 + (j-1) * Nx);
            end
            if i ~= Nx
                rhs = rhs + horizontal_term * x(i+1 + (j-1) * Nx);
            end
            if j ~= 1
                rhs = rhs + vertical_term * x(i + (j-2) * Nx);
            end
            if j ~= Ny
                rhs = rhs + vertical_term * x(i + (j) * Nx);
            end
            residual_vec(i + (j-1) * Nx) = residual_vec(i + (j-1) * Nx) + main_term * x(i + (j-1) * Nx) + rhs;
        end
    end
  
    % norm(.) gives the 2-norm of the error vector. To use the residual
    % defined in the workshet, we must divide by sqrt(N)
    residual = norm(residual_vec) / sqrt(nb_dofs);
    n = n + 1;
    
    x = reshape(x, Nx, Ny);
end

    if isinf(residual) || isnan(residual) || residual > tol_residual
        err.message = 'Could not find a solution';
        err.identifier = 'gauss_seidel_method:divergence';
        error(err);
    end

% TODO
% Implement storage double precision HERE

end