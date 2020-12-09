close all;
clear all;
clc;

% Defining the ODE and the exact solution
ode = @(p) 7 * (1 - p/10) * p;
solution = @(t) 200 ./ (20 - 10 * exp(-7 * t));

% Initial parameters
t_end = 5;
delta_t = [0.5 0.25 0.125 0.0625 0.03125];
p0 = 20;

% Cell array where results(i, j) is the array of points for the i-th
% delta_t computed with the method j (Euler = 1, Heun = 2)
results = cell(length(delta_t), 2);

% Exact error for the i-th delta_t and j-th method
error_exact = zeros(length(delta_t), 2);

for i = 1:length(delta_t)
    t = 0:delta_t(i):t_end;
    
    results{i, 1} = odeEULER(ode, p0, delta_t(i), t_end);
    results{i, 2} = odeHEUN(ode, p0, delta_t(i), t_end);
    
    % Calling exact value of the analytical solution
    p_exact = solution(t.');
    
    % Calculating the exact error for each method in the current time step
    for j = 1:2
        error_exact(i, j) = sqrt(delta_t(i) / t_end * sum((results{i, j} - p_exact).^2));
    end
    
    legend_info{i} = sprintf('Step size = %1.3f', delta_t(i));
    
    figure(1);
    plot(t, results{i, 1}(:, 1), 'LineWidth', 1);
    hold on;
    
    figure(2);
    plot(t, results{i, 2}(:, 1), 'LineWidth', 1);
    hold on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING UP THE FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titles = ["Euler method" "Heun method"];
tt = linspace(0, t_end);
analytical_solution = solution(tt);

for i = 1:2
    figure(i);

    title(titles(i));
    xlabel('Time');
    ylabel('Population');
    grid on;
    box on;
    
    xlim([0 5])
    ylim([0 20])
    
    plot(tt, analytical_solution, '--', 'LineWidth', 2);
    
    legend_info{length(delta_t) + 1} = ('Analytical solution');
    legend(legend_info);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating the factor of error reduction if the step size is halved.
% Error reduction is the error of a given time step divided by the error of
% the previous, bigger, time step, hence an element-wise division of the
% array with itself but "shifted". A padding row is added because the first
% time step does not have any previous time step to compare the error with.

% The element at error_reduction[i, j] gives the error reduction from the
% i-1 to the i step, with the j-th method.
error_reduction = [0, 0; error_exact(2:end, :) ./ error_exact(1:end-1, :)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating the approximate error
error_app = zeros(length(delta_t), 2);

sample_t = 0:delta_t(end):t_end;

for i = 1:(length(delta_t) - 1)
    t = 0:delta_t(i):t_end; 
    
    % To compute the approximated error, we need to compare solutions at
    % different time steps with our best approximation. This results in
    % comparing solutions with a different number of discrete points. In
    % order to match the number of points, we use `interp1` to make our
    % best approximation have the same number of points as the one we are
    % comparing it with. Since the time steps are factors of eachother,
    % this does not change the value of the points, but instead only
    % reduces the number of points. This is also more flexible in the case
    % where the time steps are not factors of eachother.
    
    delta_euler = results{i, 1}.' - interp1(sample_t, results{end, 1}, t);
    delta_heun = results{i, 2}.' - interp1(sample_t, results{end, 2}, t);
    
    error_app(i, 1) = sqrt(delta_t(i) / t_end * dot(delta_euler, delta_euler));
    error_app(i, 2) = sqrt(delta_t(i) / t_end * dot(delta_heun, delta_heun));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Displaying errors
T_EULER = table( ...
    [delta_t; error_exact(:, 1).'; error_reduction(:, 1).'; error_app(:, 1).'], ...
    'VariableNames', {'Explicit Euler method (q = 1)'}, ...
    'RowNames', {'δt', 'error', 'error red.', 'error app.'} ...
)

T_HEUN = table( ...
    [delta_t; error_exact(:, 2).'; error_reduction(:, 2).'; error_app(:, 2).'], ...
    'VariableNames', {'Method of Heun (q = 2)'}, ...
    'RowNames', {'δt', 'error', 'error red.', 'error app.'} ...
)

