% Runge-Kutta method (4th order) calculation
function [t,p] = odeRUNGE(ODE,t_last,delta_t,p_init)
t_init = 0;
N = (t_last-t_init)/delta_t;               % Number of calculated points
p = zeros(N,1);                            % Initializing p vector
t = zeros(N,1);                            % Initializing time vector
t(1) = 0;                                  % Initial value for time
p(1) = p_init;                             % Initial value for p
for i=1:N-1
    t(i+1)=t(i) + delta_t;
    Y1 = ODE(p(i));
    Y2 = ODE(p(i)+0.5*Y1*delta_t);
    Y3 = ODE(p(i)+0.5*Y2*delta_t);
    Y4 = ODE(p(i)+Y3*delta_t);
    
    % Definition: y(n+1)=y(n)+1/6(Y1+2(Y2)+2(Y3)+Y4)Δt
    p(i+1) = p(i)+((Y1+2*Y2+2*Y3+Y4)/6)*delta_t;
end
end