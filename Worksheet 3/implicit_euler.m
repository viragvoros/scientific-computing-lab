function p = implicit_euler(ODE, p0, delta_t, t_end, dODE)
% implicit_euler uses the Implicit Euler Method to compute discrete points of
% the solution function for the given ODE.
%
% Inputs:
%   ODE     = Derivative function that takes the p value.
%   p0      = Initial value (at t0 = 0).
%   delta_t = Size of the time step.
%   t_end   = Last t value where to stop stepping.
%   dODE    = Derivative of the ODE, for Newton's Method.
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
    % Definition: y(n+1) = y(n) + y'(n+1) * Δt
    % We change it to 0 = y(n) + y'(n+1) * Δt - y(n+1) and find the root,
    % which is the value of p(i+1). The starting value is set to p(i)
    % because it's a reasonably close value to p(i+1).
    func = @(y) p(i) + ODE(y) * delta_t - y;
    dfunc = @(y) dODE(y) * delta_t - 1;
    
    try
        p(i+1) = newton_method(func, dfunc, p(i), 1e-4, 1000);
    catch err
        fprintf('Implicit Euler divergence at step %d with dt %f.', i, delta_t);
        p = NaN;
        return;
    end
end

end