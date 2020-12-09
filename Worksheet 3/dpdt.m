function dp = dpdt(p)
% dpdt computes the derivative represented by the ODE defined as
% dp/dt = 7 * (1 - p/10) * p.
%
% Inputs:
%   p  = Value of p in the ODE.
%
% Outputs:
%   dp = Derivative value at p.

dp = 7 * (1 - p/10) * p;

end
