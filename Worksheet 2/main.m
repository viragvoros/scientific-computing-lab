
% TODO
% Calculating error approximation C)iv)
% E_approx(idelta_t)=sqrt((idelta_t/t_last)*sum((p-p_best).^2)); % error with p_best
% Probably it should be in all the three for loops


clear all
close all
clc

% Initial parameters
t_init=0;                                       % Initial value for time
t_last=5;                                       % Last value of time 
delta_t=[1 0.5 0.25 0.125];                     % Timestep sizes
p_init=1;                                       % Initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EULER method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling Euler method with different values for step size 
for idelta_t=1:length(delta_t) 
[t,p_EULER]=odeEULER(@dpdt,t_init,t_last,delta_t(idelta_t),p_init); % t,p output variables

% Calling exact value of the analytical solution
[p_exact]=calcEXACT(t);

% Calculating error C)ii)
E_EULER(idelta_t)=sqrt((idelta_t/t_last)*sum((p_EULER-p_exact).^2));

% Plotting Euler method with different values for step size
figure(1)
plot(t,p_EULER(:,1),'LineWidth',2)
xlabel('Time')
ylabel('Population')
title('Euler method')
grid on
legend_info{idelta_t}=sprintf('Step size = %1.3f',delta_t(idelta_t));
box on
hold on
end

% Plotting exact value of the analytical solution
plot(t,p_exact(:,1),'LineWidth',2)
legend_info{length(delta_t)+1}=('Analytical solution');
legend(legend_info);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEUN method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling Heun method with different values for step size 
for idelta_t=1:length(delta_t) 
[t,p_HEUN]=odeHEUN(@dpdt,t_init,t_last,delta_t(idelta_t),p_init); % t,p output variables

% Calling exact value of the analytical solution
[p_exact]=calcEXACT(t);

% Calculating error C)ii)
E_HEUN(idelta_t)=sqrt((idelta_t/t_last)*sum((p_HEUN-p_exact).^2));

% Plotting Heun method with different values for step size
figure(2)
plot(t,p_HEUN(:,1),'LineWidth',2)
xlabel('Time')
ylabel('Population')
title('Heun method')
grid on
box on
hold on
end

% Plotting exact value of the analytical solution
plot(t,p_exact(:,1),'LineWidth',2)
legend(legend_info);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUNGE-KUTTA method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling Runge-Kutta method with different values for step size 
for idelta_t=1:length(delta_t) 
[t,p_RUNGE]=odeRUNGE(@dpdt,t_init,t_last,delta_t(idelta_t),p_init); % t,p output variables

% Calling exact value of the analytical solution
[p_exact]=calcEXACT(t);

% Calculating error C)ii)
E_RUNGE(idelta_t)=sqrt((idelta_t/t_last)*sum((p_RUNGE-p_exact).^2));

% Plotting Runge-Kutta method with different values for step size
figure(3)
plot(t,p_RUNGE(:,1),'LineWidth',2)
xlabel('Time')
ylabel('Population')
title('Runge-Kutta method')
grid on
box on
hold on
end

% Plotting exact value of the analytical solution
plot(t,p_exact(:,1),'LineWidth',2)
legend(legend_info);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating factor if step size is halved C)iii)
factor_EULER = [0 E_EULER(2)/E_EULER(1) E_EULER(2)/E_EULER(3) E_EULER(3)/E_EULER(4)];
factor_HEUN = [0 E_HEUN(2)/E_HEUN(1) E_HEUN(2)/E_HEUN(3) E_HEUN(3)/E_HEUN(4)];
factor_RUNGE = [0 E_RUNGE(2)/E_RUNGE(1) E_RUNGE(2)/E_RUNGE(3) E_RUNGE(3)/E_RUNGE(4)];
%--------------------------------------------------------------------------

% Displaying errors
errors_EULER=[delta_t; E_EULER; factor_EULER]
errors_HEUN=[delta_t; E_HEUN; factor_HEUN]
errors_RUNGE=[delta_t; E_RUNGE; factor_RUNGE]
%--------------------------------------------------------------------------


