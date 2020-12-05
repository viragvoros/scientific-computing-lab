clear all;
close all;
clc;

% Initial parameters
t_end = 5;                                        % Last value of time
delta_t = [1 0.5 0.25 0.125];                     % Timestep sizes
p0 = 1;                                           % Initial condition

% Cell array where results(i, j) is the array of points for the i-th delta_t 
% computed with the method j (Euler = 1, Heun = 2, Runge-Kutta = 3)
results = cell(length(delta_t), 3);

for i = 1:length(delta_t)
    t = 0:delta_t(i):t_end;
    
    p_EULER = odeEULER(@dpdt, p0, delta_t(i), t_end);
    p_HEUN = odeHEUN(@dpdt, p0, delta_t(i), t_end);
    p_RUNGE = odeRUNGE(@dpdt, p0, delta_t(i), t_end);
    
    %% TODO : MERGE
    results{i, 1} = p_EULER;
    %% 
    
    % Calling exact value of the analytical solution
    p_exact = calcEXACT(t.');
    
    % Calculating error C)ii)
    E_EULER(i) = sqrt((delta_t(i) / t_end) * sum((p_EULER - p_exact).^2));
    E_HEUN(i) = sqrt((delta_t(i) / t_end) * sum((p_HEUN - p_exact).^2));
    E_RUNGE(i) = sqrt((delta_t(i) / t_end) * sum((p_RUNGE - p_exact).^2));
    
    legend_info{i} = sprintf('Step size = %1.3f', delta_t(i));
    
    figure(1);
    plot(t, p_EULER(:, 1), 'LineWidth', 1);
    hold on;
    
    figure(2);
    plot(t, p_HEUN(:, 1), 'LineWidth', 1);
    hold on;
    
    figure(3);
    plot(t, p_RUNGE(:, 1), 'LineWidth', 1);
    hold on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING UP THE FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titles = ["Euler method" "Heun method" "Runge-Kutta method"];
tt = linspace(0, t_end);
analytical_solution = calcEXACT(tt);

for i = 1:3
    figure(i);

    title(titles(i));
    xlabel('Time');
    ylabel('Population');
    grid on;
    box on;
    
    plot(tt, analytical_solution, '--', 'LineWidth', 2);
    
    legend_info{length(delta_t) + 1} = ('Analytical solution');
    legend(legend_info);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating factor if step size is halved C)iii)
% Error reduction is error for a given time step divided by error with the previous,
% bigger, time step, hence a component-wise division of the array with itself but "shifted"
% A padding 0 is added to have the same size as the number of time steps,
% to make it possible to write everything in a table
error_reduction_EULER = [0, E_EULER(2:end) ./ E_EULER(1:end-1)];
error_reduction_HEUN = [0, E_HEUN(2:end) ./ E_HEUN(1:end-1)];
error_reduction_RUNGE = [0, E_RUNGE(2:end) ./ E_RUNGE(1:end-1)];
%--------------------------------------------------------------------------

% Displaying errors
errors_EULER = [delta_t; E_EULER; error_reduction_EULER]
errors_HEUN = [delta_t; E_HEUN; error_reduction_HEUN]
errors_RUNGE = [delta_t; E_RUNGE; error_reduction_RUNGE]
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the approximate error C)iv)
errors_red_EULER = zeros(1, length(delta_t));

sample_t = 0:delta_t(end):t_end;

for i = 1:(length(delta_t)-1)
    t = 0:delta_t(i):t_end; 
    
    delta_euler = (results{i, 1}' - interp1(sample_t, results{length(delta_t), 1}, t));
    error_euler = sqrt (dot(delta_euler, delta_euler) * delta_t(i) / t_end)
end