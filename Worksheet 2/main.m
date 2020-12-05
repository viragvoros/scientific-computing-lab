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
%Approximation error for the i-th delta_t and j-th method
ERRORS = zeros(length(delta_t), 3);

for i = 1:length(delta_t)
    t = 0:delta_t(i):t_end;
    
    results{i, 1} = odeEULER(@dpdt, p0, delta_t(i), t_end);
    results{i, 2} = odeHEUN(@dpdt, p0, delta_t(i), t_end);
    results{i, 3} = odeRUNGE(@dpdt, p0, delta_t(i), t_end);
    
    
    % Calling exact value of the analytical solution
    p_exact = calcEXACT(t.');
    
    % Calculating error C)ii)
    
    for j=1:3
        ERRORS(i, j) = sqrt((delta_t(i) / t_end) * sum((results{i, j} - p_exact).^2));
    end
    E_EULER(i) = sqrt((delta_t(i) / t_end) * sum((results{i, 1} - p_exact).^2));
    E_HEUN(i) = sqrt((delta_t(i) / t_end) * sum((results{i, 2} - p_exact).^2));
    E_RUNGE(i) = sqrt((delta_t(i) / t_end) * sum((results{i, 3} - p_exact).^2));
    
    legend_info{i} = sprintf('Step size = %1.3f', delta_t(i));
    
    figure(1);
    plot(t, results{i, 1}(:, 1), 'LineWidth', 1);
    hold on;
    
    figure(2);
    plot(t, results{i, 2}(:, 1), 'LineWidth', 1);
    hold on;
    
    figure(3);
    plot(t, results{i, 3}(:, 1), 'LineWidth', 1);
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

% There is transposition such that error_reduction[i,j] gives the error
% reduction from the i-1 to the i step, with the j-th method
error_reduction = [[0,0,0]', ERRORS(2:end, :) ./ ERRORS(1:end-1, :)]';

error_reduction_EULER = [0, E_EULER(2:end) ./ E_EULER(1:end-1)];
error_reduction_HEUN = [0, E_HEUN(2:end) ./ E_HEUN(1:end-1)];
error_reduction_RUNGE = [0, E_RUNGE(2:end) ./ E_RUNGE(1:end-1)];

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the approximate error C)iv)
errors_app_EULER = zeros(1, length(delta_t));
errors_app_HEUN = zeros(1, length(delta_t));
errors_app_RUNGE = zeros(1, length(delta_t));


sample_t = 0:delta_t(end):t_end;

for i = 1:(length(delta_t)-1)
    t = 0:delta_t(i):t_end; 
    
    delta_euler = (results{i, 1}' - interp1(sample_t, results{length(delta_t), 1}, t));
    delta_heun = (results{i, 2}' - interp1(sample_t, results{length(delta_t), 2}, t));
    delta_runge = (results{i, 3}' - interp1(sample_t, results{length(delta_t), 3}, t));
    
    errors_app_EULER(i) = sqrt (dot(delta_euler, delta_euler) * delta_t(i) / t_end);
    errors_app_HEUN(i) = sqrt (dot(delta_heun, delta_heun) * delta_t(i) / t_end);
    errors_app_RUNGE(i) = sqrt (dot(delta_heun, delta_heun) * delta_t(i) / t_end);
end

%--------------------------------------------------------------------------

% Displaying errors
T_EULER = table([delta_t; E_EULER; error_reduction_EULER; errors_app_EULER],'VariableNames',{'explicit Euler method (q = 1)'}, 'RowNames',{'δt','error','error red.','error app.'})
T_HEUN = table([delta_t; E_HEUN; error_reduction_HEUN; errors_app_HEUN],'VariableNames',{'method of Heun (q = 2)'}, 'RowNames',{'δt','error','error red.','error app.'})
T_RUNGE = table([delta_t; E_RUNGE; error_reduction_RUNGE; errors_app_RUNGE],'VariableNames',{'Runge-Kutta method (q = 4)'}, 'RowNames',{'δt','error','error red.','error app.'})

