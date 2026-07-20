function plot_region_surface(fignum, r_grid, t_vec, C, R_cell, region_name, species_name, cmax_val)
    figure(fignum);
    surf(r_grid, t_vec, C');
    shading interp;
    caxis([0 cmax_val]);
    xlabel('r (\mum)'); ylabel('t (s)'); zlabel('C (\muM)');
    title(sprintf('%s Concentration - %s Region', species_name, region_name));
    view(45, 25);
    colorbar;
    hold on;
    xline(R_cell, 'r-', 'LineWidth', 2, 'Label', 'cell boundary');
    hold off;
end