Nx = 200;
Ny = 200;
hx = 1/(Nx+1);
hy = 1/(Ny+1);

b = zeros(Nx, Ny);
b_func = @(x,y) -2*pi^2 * sin(pi * x) * sin(pi * y);
for i = 1:Nx
   for j = 1:Ny
    b(i,j) = b_func(i*hx, j*hy);
   end
end

times = zeros(1, 3);
x = zeros(Nx+2, Ny+2);

%G-S time
tic;
x(2:end-1, 2:end-1) = gauss_seidel_solve(Nx, Ny, b);
times(1) = toc;

%Full matrix time
A = construct_matrix(Nx, Ny);
b_flat = reshape(b, Nx*Ny, 1);
tic;
x(2:end-1, 2:end-1) = reshape(A\b_flat, Nx, Ny);
times(2) = toc;

%Sparse matrix
A = sparse(A);
tic;
x(2:end-1, 2:end-1) = reshape(A\b_flat, Nx, Ny);
times(3) = toc;

fprintf("Gauss-Seidel : %f s, full matrix : %f s, sparse matrix : %f s", times(1), times(2), times(3));

% plot (contour or contourf ? COntourf seems prettier but WS says contour)
%contourf(x);