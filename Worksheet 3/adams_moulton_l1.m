function p = adams_moulton_l1(p0, delta_t, t_end)
% adams_moulton_l1 uses the Adams-Moulton Method in a linearized form to
% compute discrete points of the solution function for the given problem
% from the worksheet, even when Newton method cannot find solutions.
%
% Inputs:
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
    % Definition:
    % y(n+1) = y(n) + (7 * (1 - y(n) / 10) * y(n) + 7 * (1 - y(n+1) / 10) * y(n)) * Δt / 2
    %
    % We change it to find the value of p(i+1) in an explicit form:
    % y(n+1) = (y(n) + Δt / 2 * ((7 * (1 - y(n) / 10) * y(n)) + 7 * y(n)))...
    %          / (1 + 7 * y(n) * Δt / 20)
    
    p(i+1) = (p(i) + delta_t / 2 * ((7 * (1 - p(i) / 10) * p(i)) + 7 * p(i)))...
             / (1 + 7 * p(i) * delta_t / 20);
    
end

end