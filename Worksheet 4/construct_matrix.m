function A = construct_matrix(Nx, Ny)
% Computes the matrix for a 2D laplacian discretization of the unit square with 
% homogenous Dirichlet BC (i.e. all nodes are internal).
% Inputs :
%   Nx : number of inner nodes in the horizontal direction
%   Ny : number of inner nodes in the vertical direction
% Outputs :
%   A square matrix of dimension (Nx*Ny) where the laplacien at point (x_i,
%   y_j) is expressed on the line i + (j-1) * Nx. DOFs must be ordered
%   accordingly


hx = 1/(Nx+1);
hy = 1/(Ny+1);

nb_dofs = Nx * Ny;

A = zeros(nb_dofs);
% Dof (i, j) gets mapped to 1D index "i + (j-1) * Nx"

for i = 1:Nx
   for j = 1:Ny
       global_dof = i + (j-1) * Nx;
       A(global_dof, global_dof) = -2/hx^2 -2/hy^2;
       if i ~= 1
            neighbor_dof = (i-1) + (j-1) * Nx;
            A(neighbor_dof, global_dof) = 1/hx^2;
       end
       if i ~= Nx
            neighbor_dof = (i+1) + (j-1) * Nx;
            A(neighbor_dof, global_dof) = 1/hx^2;
       end
       if j ~= 1
            neighbor_dof = i + (j-2) * Nx;
            A(neighbor_dof, global_dof) = 1/hy^2;
       end
       if j ~= Ny
            neighbor_dof = i + (j) * Nx;
            A(neighbor_dof, global_dof) = 1/hy^2;
       end
   end
end

end