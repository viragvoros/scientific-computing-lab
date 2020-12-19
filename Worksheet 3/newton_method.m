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

if isinf(x) || isnan(x) || abs(f(x)) > tol
   err.message = 'Couldn'' find a precise solution';
   err.identifier = 'newton_method:divergence';
   error(err);
end

end