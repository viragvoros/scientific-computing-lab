function [t, p] = odeRUNGE(ODE, p0, delta_t, t_last)
% odeRUNGE uses the Fourth Order Runge-Kutta Method (RK4) to compute
% discrete points of the solution function for the given ODE.
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
    
    k1 = ODE(p(i));
    k2 = ODE(p(i) + 0.5 * k1 * delta_t);
    k3 = ODE(p(i) + 0.5 * k2 * delta_t);
    k4 = ODE(p(i) + k3 * delta_t);
    
    % Definition: y(n+1) = y(n) + (1/6)(k1 + 2(k2) + 2(k3) + k4) * Δt
    p(i+1) = p(i) + ((k1 + 2*k2 + 2*k3 + k4) / 6) * delta_t;
end

end