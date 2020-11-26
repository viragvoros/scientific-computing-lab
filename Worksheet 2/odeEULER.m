% Euler method calculation
function [t,p] = odeEULER(ODE,t_init,t_last,delta_t,p_init)
N = (t_last-t_init)/delta_t;                % Number of calculated points
p = zeros(N,1);                             % Initializing p vector
t = zeros(N,1);                             % Initializing time vector
t(1) = 0;                                   % Initial value for time
p(1) = p_init;                              % Initial value for p
for i=1:N-1
    t(i+1)=t(i) + delta_t;
    
    % Definition: y(n+1)=y(n)+y'(n)Δt)
    p(i+1)=p(i) + ODE(t(i),p(i))*delta_t;
end
end
 
