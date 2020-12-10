function x = newton_method(f, df, x0, tol)
% Newton method to solve the equation f(x) = 0
%
% Inputs:
%   f     = function to find the root of
%   df      = derivative of f
%   x0 = Initial guess
%   tol   = Accuracy condition. The iterations stop when |f(x)| < tol
%
% Outputs:
%   x       = Approximate solution

x = x0;
nb_iterations = 0;
nb_iterations_max = 100; %Arbitrary parameter

while (abs(f(x)) >= tol && nb_iterations < nb_iterations_max)
   x = x - f(x) /  df(x);
   nb_iterations = nb_iterations + 1;
end

if (nb_iterations >= nb_iterations_max || isinf(x) || isnan(x))
   error("Newton method doesn't converge ! Try another guess value or tolerance") 
end

end