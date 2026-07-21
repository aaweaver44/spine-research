function plot_snapshot_grid(fignum, C_cell, C_ecm, r_cell_grid, r_ecm_grid, ...
                            R_cell, t_vec, plot_every_unused, species_name, cmax_val)
    r_full = [0, r_cell_grid, r_ecm_grid(2:end)];
    N_theta = 200;
    theta = linspace(0, 2*pi, N_theta);
    [R_plot, Theta] = meshgrid(r_full, theta);
    X = R_plot .* cos(Theta);
    Y = R_plot .* sin(Theta);

    figure(fignum);
    set(gcf, 'Color', 'w');
    sgtitle(sprintf('%s — Cross-Section Over Time', species_name), ...
            'FontSize', 16, 'FontWeight', 'bold');

    num_snapshots = 6;
    nT = length(t_vec);
    % linear spacing
    %snapshot_indices = round(linspace(1, size(sol, 2), num_snapshots));
    % logarithmic spacing
    snapshot_indices = unique(round(logspace(0, log10(nT), num_snapshots)));

    for s = 1:length(snapshot_indices)
        j = snapshot_indices(s);
        C_full = [C_cell(1,j), C_cell(:,j)', C_ecm(2:end,j)'];
        C_mesh = repmat(C_full, N_theta, 1);
        subplot(2, 3, s);
        pcolor(X, Y, C_mesh);
        shading interp;  colormap(jet);  caxis([0 cmax_val]);
        axis equal tight off;
        hold on;
        plot(R_cell*cos(theta), R_cell*sin(theta), 'w-', 'LineWidth', 1.5);
        hold off;
        title(sprintf('t = %.3f s', t_vec(j)), 'FontSize', 12);
    end
    cb = colorbar('Position', [0.93, 0.15, 0.02, 0.7]);
    cb.Label.String = sprintf('%s conc.', species_name);
end