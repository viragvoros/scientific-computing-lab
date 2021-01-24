% Solving 2D heat equation explicitly and implicitly with Gauss-Seidel

clear all
close all
clc

L = 1;                      % domain length along x and y                                          
n = [3, 7, 15, 31];         % number of grid points along x and y
tolerance = 1e-6;           % convergence criterion
final_t = [1/8, 2/8, 3/8, 4/8]; % desired time of temperature distribution
dt = [1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]; % time step size

% Examination of stability is based on the physical meaning of the
% calculated population values. Stability denotes the applicability of a
% method. Our simple criterion was that negative values are non-applicable,
% therefore unstable. Stable cases are marked by a cross in the
% table_stable_cases tabular.

% String array stability(l, m, t) is the array of stability for the explicit
% method for the l-th δt, the m-th grid size and the t-th desired final time.
stability = string(zeros(length(dt), length(n), length(final_t)));

for t = 1:length(final_t)
for m = 1:length(n)
    for l = 1:length(dt)
        % explicit scheme
        [x_exp,y_exp,T_exp] = explicit(L, n(m), final_t(t), dt(l));
        
        %testing stability criterion        
        if all(T_exp >= 0)
          stability(l,m,t) = 'x';
        else
          stability(l,m,t) = ' ';
        end
        
        % creating a surface plot
        figure(t)
        subplot(length(n), length(dt), l+(m-1)*length(dt))
        surf(x_exp,y_exp,T_exp);
        xlabel('X axis');
        ylabel('Y axis');
        title_text_exp = sprintf(['Nx = Ny = %d \nδt = %g'], n(m), dt(l));
        title(title_text_exp);
        sgtitle_text_exp = sprintf(['Explicitly at t = %g'], final_t(t));
        sgtitle(sgtitle_text_exp);
    end
    
    % implicit gauss-seidel
    try
    [x_imp,y_imp,T_imp] = implicit_gauss_seidel(L, n(m), tolerance, final_t(t), dt(1));
    catch err
        fprintf('Gauss-Seidel divergence at Nx %d Ny %d and t = %g.', n(m), n(m), final_t(t));
        x = NaN;    
    end
    % creating a surface plot
    figure(length(final_t)+1)
    subplot(length(n), length(final_t), t+(m-1)*length(final_t))
    surf(x_imp,y_imp,T_imp);
    xlabel('X axis');
    ylabel('Y axis');
    title_text_imp = sprintf(['Nx = Ny = %d \nt = %g'], n(m), final_t(t));
    title(title_text_imp);
    sgtitle_text = sprintf(['Implicitly using Gauss-Seidel δt = %g'], dt(1));
    sgtitle(sgtitle_text);
end
end

for t = 1:length(final_t)
table_stable_cases = array2table(...
    stability(:,:,t)', ...
    'VariableNames', { ...
        'δt = 1/64', ...
        'δt = 1/128', ...
        'δt = 1/256', ...
        'δt = 1/512', ...
        'δt = 1/1024', ...
        'δt = 1/2048', ...
        'δt = 1/4096', ...
    }, ...
    'RowNames', { ...
        'Nx = Ny = 3', ...
        'Nx = Ny = 7', ...
        'Nx = Ny = 15', ...
        'Nx = Ny = 31', ...
    });

fprintf(['\n\n\nExplicitly at t = %g\n\n'], final_t(t));
disp(table_stable_cases);
end

