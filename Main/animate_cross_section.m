function animate_cross_section(fignum, C_cell, C_ecm, r_cell_grid, r_ecm_grid, ...
                               R_cell, t_vec, dt, plot_every, species_name, cmin_val, cmax_val)
 
    r_full = [0, r_cell_grid, r_ecm_grid(2:end)];
    theta  = linspace(0, 2*pi, 200);
    [R_plot, Theta] = meshgrid(r_full, theta);
    X = R_plot .* cos(Theta);
    Y = R_plot .* sin(Theta);
    x_circle = R_cell * cos(theta);
    y_circle = R_cell * sin(theta);

    figure(fignum);
    for j = 1:length(t_vec)
        if mod(j-1, plot_every) ~= 0
            continue
        end
        C_full = [C_cell(1,j), C_cell(:,j)', C_ecm(2:end,j)'];
        C_2D   = repmat(C_full, length(theta), 1);
        pcolor(X, Y, C_2D);
        shading interp;  colorbar;  colormap(jet);
        clim([cmin_val cmax_val]);           % data-driven, per species
        axis equal;
        xlabel('x (\mum)'); ylabel('y (\mum)');
        title(sprintf('%s Concentration at t = %.3f s', species_name, t_vec(j)));
        hold on;
        plot(x_circle, y_circle, 'w-', 'LineWidth', 1);
        text(0, R_cell+5, 'cell boundary', 'Color', 'w', 'HorizontalAlignment', 'center');
        hold off;
        drawnow;
        pause(0.01);
    end
end