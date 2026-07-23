function print_diagnostics(dr, dt, Nt, t_end, specs, fields)
%% Numerical health report for the two-region reaction-diffusion model.
%
%   specs  : struct array, one per species-region, with fields
%              .name  (char)  e.g. 'N cell'
%              .D     (num)   diffusivity   (um^2/s)
%              .k     (num)   reaction rate (1/s)
%              .L     (num)   region length (um)
%   fields : struct array, one per solution array, with fields
%              .name  (char)
%              .C     (matrix) solution to health-check
%
%   Flags the two Crank-Nicolson failure modes:
%     - reaction oscillation  when |k*dt| > 2
%     - diffusion ringing     when sigma is large (amp -> -1)
%   and reports whether each species is time-resolved or should be QSS.

fprintf('\n========================================================\n');
fprintf('  NUMERICAL DIAGNOSTICS\n');
fprintf('========================================================\n');
fprintf('  dr = %.3g um,  dt = %.4g s,  Nt = %d,  t_end = %.4g s\n', ...
        dr, dt, Nt, t_end);

%% ---- per species-region stability ----
fprintf('\n  %-10s %9s %9s %9s %9s  %s\n', ...
        'field','sigma','k*dt','amp_diff','amp_rxn','verdict');
fprintf('  %s\n', repmat('-',1,66));
for i = 1:numel(specs)
    s = specs(i);
    sigma = s.D*dt/dr^2;
    kdt   = s.k*dt;
    amp_d = (1 - 2*sigma)/(1 + 2*sigma);        % highest grid mode
    amp_r = (1 + kdt/2)/(1 - kdt/2);            % reaction mode

    v = '';
    if amp_r < 0,      v = [v 'RXN-OSC ']; end
    if amp_d < -0.9,   v = [v 'RINGS '];   end
    if isempty(v),     v = 'ok';           end

    fprintf('  %-10s %9.2f %9.3f %+9.4f %+9.4f  %s\n', ...
            s.name, sigma, kdt, amp_d, amp_r, v);
end
fprintf('\n  RXN-OSC : |k*dt| > 2, reaction term sign-flips each step\n');
fprintf('            -> need dt < 2/|k|, or implicit/QSS treatment\n');
fprintf('  RINGS   : sigma large, shortest wavelength barely damps\n');
fprintf('            -> transient unresolved; QSS or backward Euler\n');

%% ---- timescale check: resolved or QSS? ----
fprintf('\n  %-10s %12s %10s  %s\n','field','tau=L^2/D','dt/tau','regime');
fprintf('  %s\n', repmat('-',1,52));
for i = 1:numel(specs)
    s = specs(i);
    tau = s.L^2/s.D;
    ratio = dt/tau;
    if ratio > 3
        reg = 'dt >> tau : should be QSS';
    elseif ratio > 0.5
        reg = 'marginal';
    else
        reg = 'resolved';
    end
    fprintf('  %-10s %12.3g %10.3g  %s\n', s.name, tau, ratio, reg);
end

%% ---- penetration depth (diffusion-reaction balance) ----
fprintf('\n  %-10s %14s\n','field','lambda=sqrt(D/|k|)');
fprintf('  %s\n', repmat('-',1,30));
for i = 1:numel(specs)
    s = specs(i);
    if s.k ~= 0
        lam = sqrt(s.D/abs(s.k));
        fprintf('  %-10s %12.3g um   (%.1f x region length)\n', ...
                s.name, lam, lam/s.L);
    end
end

%% ---- solution health ----
fprintf('\n  %-12s %11s %11s  %s\n','field','min','max','health');
fprintf('  %s\n', repmat('-',1,52));
for i = 1:numel(fields)
    f = fields(i);
    mn = min(f.C(:));  mx = max(f.C(:));
    h = '';
    if any(~isfinite(f.C(:))),  h = [h 'NaN/Inf! ']; end
    if mn < -1e-10,             h = [h 'NEGATIVE! ']; end
    % crude oscillation detector: sign flips in consecutive spatial diffs
    % oscillation in SPACE: sign flips along the final radial profile
    d_space = diff(f.C(:,end));
    flips_s = sum(diff(sign(d_space)) ~= 0);
    if flips_s > numel(d_space)/3,  h = [h 'OSC-SPACE ']; end
    % oscillation in TIME: sign flips along a mid-radius time series
    d_time = diff(f.C(round(end/2), :));
    flips_t = sum(diff(sign(d_time)) ~= 0);
    if flips_t > numel(d_time)/3,   h = [h 'OSC-TIME ']; end

    if isempty(h),                  h = 'ok'; end
    fprintf('  %-12s %11.4g %11.4g  %s\n', f.name, mn, mx, h);
end
fprintf('========================================================\n\n');
end