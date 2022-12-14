close all;
clear all;
clc;

% Defining the ODE and the exact solution
ode = @(p) 7 * (1 - p/10) * p;
dode = @(p) 7 - 1.4 * p;
solution = @(t) 200 ./ (20 - 10 * exp(-7 * t));

% Initial parameters
t_end = 5;
delta_t = [0.5 0.25 0.125 0.0625 0.03125];
p0 = 20;

% Cell array where results(i, j) is the array of points for the i-th
% delta_t computed with the method j (Explicit Euler = 1, Heun = 2,
% Implicit Euler = 3, Adams-Moulton = 4, Adams-Moulton - l1 = 5,
% Adams-Moulton - l2 = 6).
results = cell(length(delta_t), 6);

% Exact error for the i-th delta_t and j-th method
error_exact = zeros(length(delta_t), 6);

for i = 1:length(delta_t)
    t = 0:delta_t(i):t_end;

    % For both implicit schemes (3 : implicit Euler, 4 : Adams-Moulton),
    % the solver returns NaN if no solution was found.
    results{i, 1} = explicit_euler(ode, p0, delta_t(i), t_end);
    results{i, 2} = heun(ode, p0, delta_t(i), t_end);
    results{i, 3} = implicit_euler(ode, p0, delta_t(i), t_end, dode);
    results{i, 4} = adams_moulton(ode, p0, delta_t(i), t_end, dode);
    results{i, 5} = adams_moulton_l1(p0, delta_t(i), t_end);
    results{i, 6} = adams_moulton_l2(p0, delta_t(i), t_end);
    
    % Calling exact value of the analytical solution
    p_exact = solution(t.');
    
    % Calculating the exact error for each method in the current time step
    for j = 1:6
        error_exact(i, j) = sqrt(delta_t(i) / t_end * sum((results{i, j} - p_exact).^2));
    end
    
    legend_string = sprintf('Step size = %1.3f', delta_t(i));
    
    figure(1);
    plot(t, results{i, 1}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
    hold on;
    
    figure(2);
    plot(t, results{i, 2}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
    hold on;
    
    if ~isnan(results{i, 3})
        figure(3);
        plot(t, results{i, 3}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
        hold on;
    end
    
    if ~isnan(results{i, 4})
        figure(4);
        plot(t, results{i, 4}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
        hold on;
    end
        
    figure(5);
    plot(t, results{i, 5}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
    hold on;
    
    figure(6);
    plot(t, results{i, 6}(:, 1), 'LineWidth', 1, 'DisplayName', legend_string);
    hold on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING UP THE FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titles = [ ...
    "Explicit Euler Method" ...
    "Heun Method" ...
    "Implicit Euler Method" ...
    "Adams-Moulton Method" ...
    "Adams-Moulton Method - Linearisation 1" ...
    "Adams-Moulton Method - Linearisation 2" ...
];
tt = linspace(0, t_end);
analytical_solution = solution(tt);

for i = 1:6
    figure(i);

    title(titles(i));
    xlabel('Time');
    ylabel('Population');
    grid on;
    box on;
    
    xlim([0 5])
    ylim([0 20])
    
    plot(tt, analytical_solution, '--', 'LineWidth', 2,...
        'DisplayName', 'Analytical solution');
    
    legend('Location','northeastoutside');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating the factor of error reduction if the step size is halved.
% Error reduction is the error of the previous and bigger time step divided
% by the given time step, hence an element-wise division of the array with
% itself but "shifted". A padding row is added because the first time step
% does not have any previous time step to compare the error with.

% The element at error_reduction[i, j] gives the error reduction from the
% i-1 to the i step, with the j-th method.
error_reduction = [ ...
    0, 0, 0, 0, 0, 0; ...
	error_exact(1:end-1, :) ./ error_exact(2:end, :) ...
];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Displaying errors
table_explicit_euler = table( ...
    [delta_t; error_exact(:, 1).'; error_reduction(:, 1).'], ...
    'VariableNames', {'Explicit Euler method'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

table_heun = table( ...
    [delta_t; error_exact(:, 2).'; error_reduction(:, 2).'], ...
    'VariableNames', {'Method of Heun'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

table_implicit_euler = table( ...
    [delta_t; error_exact(:, 3).'; error_reduction(:, 3).'], ...
    'VariableNames', {'Implicit Euler method'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

table_adams_moulton = table( ...
    [delta_t; error_exact(:, 4).'; error_reduction(:, 4).'], ...
    'VariableNames', {'Adams-Moulton'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

table_adams_moulton_lin1 = table( ...
    [delta_t; error_exact(:, 5).'; error_reduction(:, 5).'], ...
    'VariableNames', {'Adams-Moulton ??? linearisation 1'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

table_adams_moulton_lin2 = table( ...
    [delta_t; error_exact(:, 6).'; error_reduction(:, 6).'], ...
    'VariableNames', {'Adams-Moulton ??? linearisation 2'}, ...
    'RowNames', {'??t', 'error', 'error red.'} ...
)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Examination of stability is based on the physical meaning of the
% calculated population values. Stability denotes the applicability of a
% method. Our simple criterion was that negative population values are
% unphysical and non-applicable, therefore unstable. Stable cases are
% marked by a cross in the table_stable_cases tabular.

% String array stability(i, j) is the array of stability for the i-th
% delta_t computed with the method j (Explicit Euler = 1, Heun = 2,
% Implicit Euler = 3, Adams-Moulton = 4, Adams-Moulton - l1 = 5,
% Adams-Moulton - l2 = 6).
stability = string(zeros(length(delta_t), 6));

for i = 1:length(delta_t)
 for j = 1:6  
   if all(results{i, j} >= 0)
     stability(i,j) = 'x';
   else
     stability(i,j) = ' ';
   end
 end
end

table_stable_cases = array2table(...
    stability, ...
    'VariableNames', { ...
        'Explicit Euler', ...
        'Heun', ...
        'Implicit Euler', ...
        'Adams-Moulton', ...
        'Adams-Moulton l1', ...
        'Adams-Moulton l2' ...
     }, ...
    'RowNames', { ...
        '??t = 0.5', ...
        '??t = 0.25', ...
        '??t = 0.125', ...
        '??t = 0.0625', ...
        '??t = 0.03125' ...
    } ...
) 

