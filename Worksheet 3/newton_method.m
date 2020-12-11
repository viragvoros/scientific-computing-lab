function x = newton_method(f, df, x0, tol, max_iter)
% Newton's Method to solve the equation f(x) = 0.
%
% Inputs:
%   f        = Function to find the root of.
%   df       = Derivative of function f.
%   x0       = Initial guess.
%   tol      = Accuracy condition, the iterations stop when |f(x)| < tol.
%   max_iter = Maximum number of iterations to perform.
% 
% Outputs:
%   x        = Approximate solution.

x = x0;
nb_iterations = 0;

while abs(f(x)) >= tol && nb_iterations < max_iter
   x = x - f(x) / df(x);
   nb_iterations = nb_iterations + 1;
end

% TODO: Add back condition `nb_iterations >= max_iter` 
if isinf(x) || isnan(x)
   error("Newton method doesn't converge! Try another guess value or tolerance.") 
end

end