% Heun method calculation
function [t,p] = odeHEUN(ODE,t_last,delta_t,p_init)
t_init = 0;
N = (t_last-t_init)/delta_t;                % Number of calculated points
p = zeros(N,1);                             % Initializing p vector
t = zeros(N,1);                             % Initializing time vector
t(1) = 0;                                   % Initial value for time
p(1) = p_init;                              % Initial value for p
for i=1:N-1
    t(i+1)=t(i) + delta_t;
    Y1 = ODE(p(i));                   % y'(n)
    Y2 = ODE(p(i)+Y1*delta_t);% y'(n+1)
    
    % Definition: y(n+1)=y(n)+0.5(y'(n)+y'(n+1))Δt
    p(i+1)=p(i) + (delta_t/2)*(Y1+Y2);     
end
end