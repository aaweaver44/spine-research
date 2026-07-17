function test_cartesian_prototype
% TEST_CARTESIAN_PROTOTYPE  Verification suite for the Cartesian baseline.
%
% As of July 17 2026 THIS DOCUMENT IS 100% CLAUDE GENERATED AND UNREVIEWED
%
%   Verifies crank_nicholson_textbook.m and forward_difference_textbook.m
%   against the exact solution of the homogeneous heat equation on a slab.
%
%   Companion to: validation_01_cartesian_prototype.md
%
%   Usage:
%       cd('.../Research/Cartesian prototype')
%       test_cartesian_prototype
%
%   Prints a PASS/FAIL table. Expect ~1 minute runtime.

clc;
fprintf('\n=========================================================\n');
fprintf('  VALIDATION 01 - CARTESIAN PROTOTYPE\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('=========================================================\n');

% The solvers call surf() internally; suppress the figures they open.
old_vis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
cleanup = onCleanup(@() restore_figs(old_vis));

results = {};   % {name, pass, detail}

% Fixed problem definition
D  = 1;
L  = 1;
IC = @(x) sin(2*pi*x).^2;
BCz = @(t) 0*ones(size(t));       % constant zero (NOT the driver's lbc*t form)

%% -------------------------------------------------------------------
%  TEST 1 - Which equation does the code actually solve?
%  -------------------------------------------------------------------
fprintf('\n[Test 1] Governing equation identification\n');
M = 200; T = 0.05; N = 4000;
w  = crank_nicholson_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
x  = linspace(0, L, M+1)';
c_num  = w(:, end);
c_pure = exact_heat(x, T, D);
c_rxn  = exact_heat(x, T, D) * exp(1*T);      % if reaction k=1 were present

err_pure = max(abs(c_num - c_pure));
err_rxn  = max(abs(c_num - c_rxn));
fprintf('   max|numeric - pure diffusion|     = %.3e\n', err_pure);
fprintf('   max|numeric - diffusion+reaction| = %.3e\n', err_rxn);
p = err_pure < 1e-4 && err_rxn > 1e-3;
results(end+1,:) = {'1. Equation identification', p, ...
    sprintf('pure diffusion (err %.1e vs %.1e)', err_pure, err_rxn)};
report(p);

%% -------------------------------------------------------------------
%  TEST 2 - Spatial order of accuracy (expect 2)
%  -------------------------------------------------------------------
fprintf('\n[Test 2] Spatial order of accuracy\n');
T = 0.05; N = 20000;         % dt tiny => temporal error negligible
Ms = [10 20 40 80 160];
errs = zeros(size(Ms));
for i = 1:numel(Ms)
    M = Ms(i);
    w = crank_nicholson_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
    x = linspace(0, L, M+1)';
    errs(i) = max(abs(w(:,end) - exact_heat(x, T, D)));
end
ords = log2(errs(1:end-1) ./ errs(2:end));
for i = 1:numel(Ms)
    if i == 1
        fprintf('   M=%4d  dx=%.4f  err=%.4e\n', Ms(i), L/Ms(i), errs(i));
    else
        fprintf('   M=%4d  dx=%.4f  err=%.4e   order=%.2f\n', ...
                Ms(i), L/Ms(i), errs(i), ords(i-1));
    end
end
p = abs(ords(end) - 2) < 0.1;
results(end+1,:) = {'2. Spatial order', p, sprintf('order=%.2f (expect 2)', ords(end))};
report(p);

%% -------------------------------------------------------------------
%  TEST 3 - Temporal order, Crank-Nicolson (expect 2)
%  -------------------------------------------------------------------
fprintf('\n[Test 3] Temporal order of accuracy, Crank-Nicolson\n');
M = 1600; T = 0.05;          % M large => spatial error ~1e-7, dt dominates
Ns = [5 10 20 40 80];
errs = zeros(size(Ns));
x = linspace(0, L, M+1)';
ref = exact_heat(x, T, D);
for i = 1:numel(Ns)
    w = crank_nicholson_textbook(0, L, 0, T, M, Ns(i), D, IC, BCz, BCz, 1);
    errs(i) = max(abs(w(:,end) - ref));
end
ords = log2(errs(1:end-1) ./ errs(2:end));
for i = 1:numel(Ns)
    if i == 1
        fprintf('   N=%3d  dt=%.5f  err=%.4e\n', Ns(i), T/Ns(i), errs(i));
    else
        fprintf('   N=%3d  dt=%.5f  err=%.4e   order=%.2f\n', ...
                Ns(i), T/Ns(i), errs(i), ords(i-1));
    end
end
p = abs(ords(end) - 2) < 0.15;
results(end+1,:) = {'3. Temporal order (CN)', p, sprintf('order=%.2f (expect 2)', ords(end))};
report(p);

%% -------------------------------------------------------------------
%  TEST 4 - Forward Euler convergence at fixed sigma (expect 2)
%  -------------------------------------------------------------------
fprintf('\n[Test 4] Forward Euler convergence at fixed sigma=0.4\n');
T = 0.05; sig = 0.4;
Ms = [10 20 40 80];
errs = zeros(size(Ms));
for i = 1:numel(Ms)
    M  = Ms(i);
    dx = L/M;
    N  = round(T / (sig*dx^2/D));
    w  = forward_difference_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
    x  = linspace(0, L, M+1)';
    errs(i) = max(abs(w(:,end) - exact_heat(x, T, D)));
    fprintf('   M=%3d  dx=%.4f  N=%6d  err=%.4e\n', M, dx, N, errs(i));
end
ords = log2(errs(1:end-1) ./ errs(2:end));
fprintf('   observed orders: %s\n', num2str(ords', '%.2f  '));
p = abs(ords(end) - 2) < 0.15;
results(end+1,:) = {'4. Fixed-sigma order (FE)', p, sprintf('order=%.2f (expect 2)', ords(end))};
report(p);

%% -------------------------------------------------------------------
%  TEST 5 - Stability boundary (theory: FE stable iff sigma <= 0.5)
%  -------------------------------------------------------------------
fprintf('\n[Test 5] Stability boundary\n');
M = 10; T = 1; dts = [0.0030 0.0045 0.0050 0.0055 0.0060];
ok = true;
fprintf('     dt      sigma     FE max|c|      CN max|c|\n');
for i = 1:numel(dts)
    dt = dts(i); N = round(T/dt);
    sg = D*dt/(L/M)^2;
    wf = forward_difference_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
    wc = crank_nicholson_textbook   (0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
    mf = max(abs(wf(:))); mc = max(abs(wc(:)));
    fprintf('   %.4f   %.3f    %.4e    %.4e\n', dt, sg, mf, mc);
    % FE must be bounded below sigma=0.5 and divergent above
    if sg <= 0.5,  ok = ok && (mf < 2);
    else,          ok = ok && (mf > 10);
    end
    ok = ok && (mc < 2);          % CN bounded throughout
end
results(end+1,:) = {'5. Stability boundary', ok, 'FE diverges above sigma=0.5; CN stable'};
report(ok);

%% -------------------------------------------------------------------
%  TEST 6 - Decay rate and positivity (physical check)
%  -------------------------------------------------------------------
fprintf('\n[Test 6] Asymptotic decay rate and positivity\n');
M = 100; N = 2000; T = 1;
w  = crank_nicholson_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
mx = max(abs(w), [], 1);
tt = linspace(0, T, N+1);
monotone = all(diff(mx) <= 1e-12);
i1 = round(0.3*N); i2 = round(0.6*N);
rate = log(mx(i1)/mx(i2)) / (tt(i2) - tt(i1));
fprintf('   monotone non-increasing : %d\n', monotone);
fprintf('   measured decay rate     : %.4f\n', rate);
fprintf('   theoretical (D*pi^2)    : %.4f\n', D*pi^2);
fprintf('   relative error          : %.3f %%\n', 100*abs(rate - D*pi^2)/(D*pi^2));
fprintf('   min value (positivity)  : %.3e\n', min(w(:)));
p = monotone && abs(rate - D*pi^2)/(D*pi^2) < 0.01 && min(w(:)) > -1e-10;
results(end+1,:) = {'6. Decay rate / positivity', p, sprintf('rate=%.4f vs %.4f', rate, D*pi^2)};
report(p);

%% -------------------------------------------------------------------
%  TEST 7 - Scheme-to-scheme agreement inside the stable range
%  -------------------------------------------------------------------
fprintf('\n[Test 7] Crank-Nicolson vs Forward Euler agreement (sigma=0.3)\n');
M = 10; dt = 0.003; T = 1; N = round(T/dt);
wc = crank_nicholson_textbook   (0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
wf = forward_difference_textbook(0, L, 0, T, M, N, D, IC, BCz, BCz, 1);
dmax = max(abs(wc(:,end) - wf(:,end)));
fprintf('   max|CN - FE| at t=1 : %.3e\n', dmax);
p = dmax < 1e-4;
results(end+1,:) = {'7. Scheme agreement', p, sprintf('max diff %.1e', dmax)};
report(p);

%% -------------------------------------------------------------------
%  FINDING CHECKS - defects recorded in the report, not pass/fail
%  -------------------------------------------------------------------
fprintf('\n---------------------------------------------------------\n');
fprintf('  RECORDED FINDINGS (informational)\n');
fprintf('---------------------------------------------------------\n');

sg_default = D*0.0055/(0.1)^2;
fprintf('\n[Finding 1] Shipped default dt=0.0055, dx=0.1 -> sigma = %.3f\n', sg_default);
if sg_default > 0.5
    fprintf('   >> EXCEEDS the Forward Euler limit of 0.5.\n');
    fprintf('   >> reaction_diffusion.m as committed produces a divergent FE plot.\n');
end

fprintf('\n[Finding 2] Driver BCs are ramps: BC_L = @(t) lbc*t\n');
BCramp = @(t) 1*t;
wr = crank_nicholson_textbook(0, L, 0, 1, 10, 200, D, IC, BCramp, BCramp, 1);
fprintf('   with lbc=1:  c(0,t=0)=%.3f  c(0,t=0.5)=%.3f  c(0,t=1)=%.3f\n', ...
        wr(1,1), wr(1,101), wr(1,end));
fprintf('   >> a constant BC would hold %.3f at all t. This is a ramp.\n', 1);
fprintf('   >> harmless at lbc=0 as shipped; a trap for any nonzero value.\n');

%% -------------------------------------------------------------------
%  SUMMARY
%  -------------------------------------------------------------------
fprintf('\n=========================================================\n');
fprintf('  SUMMARY\n');
fprintf('=========================================================\n');
n_fail = 0;
for i = 1:size(results,1)
    if results{i,2}, tag = 'PASS'; else, tag = 'FAIL'; n_fail = n_fail + 1; end
    fprintf('  [%s]  %-28s  %s\n', tag, results{i,1}, results{i,3});
end
fprintf('---------------------------------------------------------\n');
fprintf('  %d/%d passed\n', size(results,1)-n_fail, size(results,1));
fprintf('=========================================================\n\n');

close all;
end

% ====================================================================
%  LOCAL FUNCTIONS
% ====================================================================

function u = exact_heat(x, t, D)
% Exact solution of u_t = D u_xx on [0,1], u(0)=u(1)=0, u(x,0)=sin^2(2 pi x).
%   u(x,t) = sum_{n odd} b_n sin(n pi x) exp(-D n^2 pi^2 t)
%   b_n    = -32 / (pi n (n^2 - 16))
% Only odd modes are excited, so n=4 (which would divide by zero) never occurs.
u = zeros(size(x));
for n = 1:2:201
    bn = -32 / (pi * n * (n^2 - 16));
    u  = u + bn * sin(n*pi*x) * exp(-D * (n*pi)^2 * t);
end
end

function report(p)
if p, fprintf('   -> PASS\n'); else, fprintf('   -> FAIL\n'); end
end

function restore_figs(v)
set(0, 'DefaultFigureVisible', v);
end
