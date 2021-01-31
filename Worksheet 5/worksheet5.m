% Solving 2D heat equation explicitly and implicitly with Gauss-Seidel

clear all
close all
clc

% Flags to enable or disable the displaying of graphs, and the saving of
% images.
display_graphs = true;
save_images = true;

L = 1;                          % domain length along x and y
n = [3, 7, 15, 31];             % number of grid points along x and y
tolerance = 1e-6;               % convergence criterion
final_t = [1/8, 2/8, 3/8, 4/8]; % desired time of temperature distribution
dt = [1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]; % time step size

T_states = cell(length(n), length(dt)+1, length(final_t), 3);

% Examination of stability is based on the physical meaning of the
% calculated population values. Stability denotes the applicability of a
% method. Our simple criterion was that negative values and values bigger
% than 1 are unphysical non-applicable, therefore unstable. Stable cases
% are marked by a cross in the tabulars of stable cases.

% String array stability(l, m, t) is the array of stability for the explicit
% method for the l-th δt, the m-th grid size and the t-th desired final time.
stability = string(zeros(length(dt), length(n), length(final_t)));

for m = 1:length(n)
    for l = 1:length(dt)
        for t = 1:length(final_t)
            
            % explicit scheme
            if t == 1
                [x_exp, y_exp, T_exp] = explicit(L, n(m), 0, final_t(t), dt(l));
            else
                [x_exp, y_exp, T_exp] = explicit(L, n(m), final_t(t-1), final_t(t), dt(l), T_exp);
            end
            
            T_states{m, l, t, 1} = x_exp;
            T_states{m, l, t, 2} = y_exp;
            T_states{m, l, t, 3} = T_exp;
            
            % testing stability criterion
            if all(T_exp(:) >= 0) && all(T_exp(:) <= 1)
                stability(l,m,t) = 'x';
            else
                stability(l,m,t) = ' ';
            end
            
        end
    end
end

for m = 1:length(n)
    for t = 1:length(final_t)
        
        % implicit gauss-seidel
        try
            if t == 1
                [x_imp, y_imp, T_imp] = implicit_gauss_seidel(L, n(m), ...
                    tolerance, 0, final_t(t), dt(1));
            else
                [x_imp, y_imp, T_imp] = implicit_gauss_seidel(L, n(m), ...
                    tolerance, final_t(t-1), final_t(t), dt(1), T_imp);
            end
        catch err
            fprintf('Gauss-Seidel divergence at Nx %d Ny %d and t = %g.', ...
                n(m), n(m), final_t(t));
            x = NaN;
        end
        
        T_states{m, length(dt)+1, t, 1} = x_imp;
        T_states{m, length(dt)+1, t, 2} = y_imp;
        T_states{m, length(dt)+1, t, 3} = T_imp;
        
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
    
    fprintf(['\n\n\nStable cases explicitly at t = %g\n\n'], final_t(t));
    disp(table_stable_cases);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting and display the graphs. This is enabled by the `display_graphs`
% flag.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if display_graphs
    for t = 1:length(final_t)
        figure(t);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        sgtitle_text_exp = sprintf(['Explicitly at t = %g'], final_t(t));
        sgtitle(sgtitle_text_exp);
        
        for m = 1:length(n)
            for l = 1:length(dt)
                subplot(length(n), length(dt), l+(m-1)*length(dt))
                surf(T_states{m, l, t, 1}, T_states{m, l, t, 2}, T_states{m, l, t, 3});
                xlabel('X');
                ylabel('Y');
                zlabel('T');
                title_text_exp = sprintf(['Nx = Ny = %d \nδt = %g'], n(m), dt(l));
                title(title_text_exp);
            end
        end
    end
    
    figure(length(final_t)+1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    sgtitle_text = sprintf(['Implicitly using Gauss-Seidel δt = %g'], dt(1));
    sgtitle(sgtitle_text);
    
    for m = 1:length(n)
        for t = 1:length(final_t)
            subplot(length(n), length(final_t), t+(m-1)*length(final_t))
            surf(T_states{m, length(dt)+1, t, 1}, ...
                T_states{m, length(dt)+1, t, 2}, T_states{m, length(dt)+1, t, 3});
            xlabel('X');
            ylabel('Y');
            zlabel('T');
            title_text_imp = sprintf(['Nx = Ny = %d \nt = %g'], n(m), final_t(t));
            title(title_text_imp);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving images of each plot in the output/ directory. This is enabled by
% the `save_images` flag. This will add around a minute or two to the
% execution time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if save_images
    mkdir('output');
    saving_figure = figure('visible', 'off');
    

    for t = 1:length(final_t)
        for m = 1:length(n)
            for l = 1:length(dt)
                set(0, 'CurrentFigure', saving_figure);
                surf(T_states{m, l, t, 1}, T_states{m, l, t, 2}, T_states{m, l, t, 3});
                xlabel('X');
                ylabel('Y');
                zlabel('T');
                sgtitle_text_exp = sprintf(['Explicitly at t = %g\nNx = Ny = %d \nδt = %g'], ...
                    final_t(t), n(m), dt(l));
                sgtitle(sgtitle_text_exp);
                exportgraphics(saving_figure, sprintf(['output/mesh_%d_dt_%g_final_t_%g.png'], ...
                    n(m), dt(l), final_t(t)));
            end
        end
    end
end
