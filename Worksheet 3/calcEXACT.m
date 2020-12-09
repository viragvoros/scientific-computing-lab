function value = calcEXACT(t)
% calcEXACT Calculates the exact value of the analytical solution defined
% as p(t) = 200 ./ (20 + 10*exp(-7 * t)).
%
% Inputs:
%   t     = Value of the variable t in the solution.
%
% Outputs:
%   value = Output of the analytical function at time t.

value = 200 ./ (20 + 10*exp(-7 * t));

end