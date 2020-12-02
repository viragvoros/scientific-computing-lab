function value = calcEXACT(t)
% calcEXACT Calculates the exact value of the analytical solution defined
% as p(t) = 10 / (1 + 9 * exp(-t)).
%
% Inputs:
%   t     = Value of the variable t in the solution.
%
% Outputs:
%   value = Output of the analytical function at time t.

value = 10 ./ (1 + 9*exp(-t));

end