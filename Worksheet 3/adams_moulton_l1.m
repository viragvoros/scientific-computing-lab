function p = adams_moulton_l1(ODE, p0, delta_t, t_end)
% adams_moulton_l1 uses the Adams-Moulton method in a linearized form to
% compute discrete points of the solution function for the given ODE, even
% for unstable time step sizes.
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
    % Definition: y(n+1) = y(n) + (y'(n) + 7 * (1 - y(n+1) / 10) * y(n)) * Δt / 2
    % We change it to 0 = y(n) + (y'(n) + 7 * (1 - y(n+1) / 10) * y(n)) * Δt / 2 - y(n+1)
    % and find the root, which is the value of p(i+1). The starting value is
    % set to p(i) because it's a reasonably close value to p(i+1).
    func = @(y) p(i) + (ODE(p(i)) + 7 * (1 - y / 10) * p(i)) * delta_t / 2 - y;
    dfunc = @(y) -7 * p(i) * delta_t / 20 - 1;
    
    p(i+1) = newton_method(func, dfunc, p(i), 1e-4, 1000);
end

end