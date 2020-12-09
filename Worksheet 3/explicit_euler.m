function p = explicit_euler(ODE, p0, delta_t, t_end)
% explicit_euler uses the Explicit Euler Method to compute discrete points of
% the solution function for the given ODE.
%
% Inputs:
%   ODE     = Derivative function that takes the p value.
%   p0      = Initial value (at t0 = 0).
%   delta_t = Size of the time step.
%   t_end   = Last t value where to stop stepping.
%
% Outputs:
%   p       = Vector containing the approximated values of p.

% Where to start stepping from
t_init = 0;

% Number of calculated points
N = floor((t_end - t_init) / delta_t) + 1;

% Initialize the p vector with the initial condition
p = zeros(N, 1);
p(1) = p0;

for i = 1:N - 1
    % Definition: y(n+1) = y(n) + y'(n) * Δt
    p(i+1) = p(i) + ODE(p(i)) * delta_t;
end

end