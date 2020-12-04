function [t, p] = odeEULER(ODE, p0, delta_t, t_last)
% odeEULER uses the Explicit Euler Method to compute discrete points of
% the solution function for the given ODE.
%
% Inputs:
%   ODE     = Derivative function that takes the p value.
%   p0      = Initial value (at t0 = 0).
%   delta_t = Size of the time step.
%   t_last  = Last t value where to stop stepping.
%
% Outputs:
%   [t, p]  = Vector containing the approximated values of p for a given t.

% Where to start stepping from
t_init = 0;

% Number of calculated points
N = floor((t_last - t_init) / delta_t) + 1;

% Initialize t and p vectors with the initial condition
t = zeros(N, 1);
p = zeros(N, 1);
t(1) = 0;
p(1) = p0;

for i = 1:N - 1
    t(i+1) = t(i) + delta_t;
    
    % Definition: y(n+1) = y(n) + y'(n) * Δt
    p(i+1) = p(i) + ODE(p(i)) * delta_t;
end

end