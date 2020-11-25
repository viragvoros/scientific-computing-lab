% EULER method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%--------------------------------------------------------------------------

% Initial parameters
t_init=0;                                       % Initial value for time
t_last=5;                                       % Last value of time 
delta_t=[1 0.5 0.25 0.125];                     % Timestep sizes
p_init=1;                                       % Initial condition
%--------------------------------------------------------------------------

% Calling Euler method with different values for step size 
for idelta_t=1:length(delta_t) 
[t,p]=odeEULER(@dpdt,t_init,t_last,delta_t(idelta_t),p_init); % t,p output variables  

% Calling analyical solution
[p_exact]=calcEXACT(t);

% Calculating error C)ii)
E(idelta_t)=sqrt((idelta_t/t_last)*sum((p-p_exact).^2));

% TO BE INTERPRETED
% Calculating error approximation C)iv)
% E_approx(idelta_t)=sqrt((idelta_t/t_last)*sum((p-p_best).^2)); % error with p_best

% Plotting Euler method with different values for step size
plot(t,p(:,1),'LineWidth',2)
xlabel('Time')
ylabel('Population')
title('Euler method')
grid on
legend_info{idelta_t}=sprintf('Step size = %1.3f',delta_t(idelta_t));
legend(legend_info);
box on
hold on
end
%--------------------------------------------------------------------------

% Calculating factor if step size is halved C)iii)
factor = [0 E(2)/E(1) E(2)/E(3) E(3)/E(4)];
%--------------------------------------------------------------------------

% Displaying error for Euler method
errors_EULER=[delta_t; E; factor]
%--------------------------------------------------------------------------

% Plotting exact value of the analytical solution
plot(t,p_exact(:,1),'LineWidth',2)
legend_info{length(delta_t)+1}=('Analytical solution');
legend(legend_info);
%--------------------------------------------------------------------------

% Calculating exact value of the analytical solution
function exact_value=calcEXACT(t)
exact_value=10./(1+9*exp(-t));
end
%--------------------------------------------------------------------------

% Defining ODE
function dpdx=dpdt(t,p)                     % ODE defined as the following 
dpdx=((1-(p/10))*p);
end 
%--------------------------------------------------------------------------

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
 
