function plot_region_surface(fignum, r_grid, t_vec, C, R_cell, region_name, species_name, cmin_val, cmax_val)
    figure(fignum);
    surf(r_grid, t_vec, C');
    shading interp;
    caxis([cmin_val cmax_val]);
    xlabel('r (\mum)'); ylabel('t (s)'); zlabel('C (\muM)');
    title(sprintf('%s Concentration - %s Region', species_name, region_name));
    view(45, 25);
    colorbar;
    hold on;
    hold off;
end