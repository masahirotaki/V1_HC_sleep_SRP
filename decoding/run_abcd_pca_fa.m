function results = run_abcd_pca_fa(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_PCA_FA  PCA and factor analysis (FA) trajectories of CA1 / V1 population
% activity during the GAVNIK ABCD sequence: single-trial and trial-averaged.
%
%   results = run_abcd_pca_fa('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Follows V1_MEC_acute_MAT/dimension_reduction/FA_continuous_test.m (DataHigh FA,
% copied to external/DataHigh):
%   - per sequence (trial): spike counts in 20 ms bins, Gaussian-smoothed (sigma 40 ms),
%     -0.5 to +1.1 s around A onset -> rates (Hz), neuron x time
%   - soft normalisation: rate / (range over all trials + 5 Hz)
%   - FA dimensionality q: 5-fold cross-validation over trials (q = 5..min(n-5, 50)),
%     q with the highest held-out log-likelihood (fastfa / fastfa_estep;
%     leave-neuron-out MSE from cosmoother_fa is stored too)
%   - FA: reducedims_EH(D, 3, q) on all trials -> orthonormalised factors (orthogonalize)
%   - PCA: reducedims_EH(D, 1, n_pc)
% The mean trajectory is the average of the single-trial projections (= projection of
% the trial-averaged activity, both maps being linear). Trial-to-trial variability:
% distance of each trial to the leave-one-out mean trajectory (top 3 dimensions).
% Model parameters (FA L / Ph / d and orthogonalising transform, PCA coeff / mean,
% softnorm scale, cluster ids) are saved for projecting sleep activity.
%
% Outputs (decoding_* names) go to <analysis>/<subject>/<date>/<stimulus>/Decoding.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
opts = fill_defaults(opts);
set(groot, 'defaultFigureVisible', cfg.plot.visible);

paths  = get_session_paths(subject, date, stimulus_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
trials = get_abcd_trials(paths);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
check_output_path(out_dir, cfg);
n_seq = numel(trials.seq_onset);
if isempty(opts.window)
    opts.window = [-0.5, round(trials.seq_duration, 2) + 0.5];
end
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));

% binning with padding so smoothing has no edge effects, then trim
pad = 4 * opts.sigma;
edges = (opts.window(1) - pad):opts.bin:(opts.window(2) + pad + 1e-9);
centres = edges(1:end - 1) + opts.bin / 2;
keep_t = centres > opts.window(1) & centres < opts.window(2);
t = centres(keep_t);
phase = phase_of(t, trials);

% DataHigh fastfa calls randn('state', 0), which puts MATLAB in legacy RNG mode, and every later
% rng(seed) in the session would then fail; restore the default generator on exit
restore_rng = onCleanup(@() rng('default'));
base = sprintf('decoding_%s', stimulus_name);
figs = gobjects(0);
names = {};
for p = 1:numel(opts.populations)
    P = opts.populations(p);
    sel = units.region == P.region & (units.is_good | ~P.good_only);
    ids = units.cluster_id(sel);
    counts = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset, edges);
    rates = smooth_counts(counts, opts) / opts.bin;          % trial x neuron x time (Hz)
    rates = rates(:, :, keep_t);
    mean_rate = squeeze(mean(rates, [1 3]));
    active = mean_rate(:) >= opts.min_rate;
    ids = ids(active);
    rates = rates(:, active, :);
    n = numel(ids);
    fprintf('\n=== %s: %d units (>= %.1f Hz in the trial windows)\n', P.name, n, opts.min_rate);
    if n < 10
        fprintf('    too few units for FA cross-validation (q from 5 to n - 5) - skipped\n');
        continue
    end

    % soft normalisation over all trial windows
    flat = reshape(permute(rates, [2 3 1]), n, []);
    scale = max(flat, [], 2) - min(flat, [], 2) + opts.softnorm_const;
    D = struct('data', {}, 'condition', {}, 'type', {});
    for s = 1:n_seq
        D(s).data = squeeze(rates(s, :, :)) ./ scale;   % neuron x time
        D(s).condition = 'ABCD';
        D(s).type = 'traj';
    end

    pop = struct('name', P.name, 'region', P.region, 'good_only', P.good_only, 'cluster_id', ids, ...
        'softnorm_scale', scale, 'mean_rate_hz', mean_rate(active)');

    %% FA: choose q by cross-validation over trials
    cv = fa_cross_validate(D, opts, n);
    q = cv.q_opt;
    fprintf('FA q_opt = %d (candidates %d-%d)\n', q, cv.candidate_dims(1), cv.candidate_dims(end));
    [newD, C, ~, explained, params] = reducedims_EH(D, 3, q, []);
    [Zall, ~] = fastfa_estep([D.data], params);
    [~, Lorth, TT] = orthogonalize(Zall.mean, params.L);
    shared = sum(params.L.^2, 2) ./ (sum(params.L.^2, 2) + params.Ph);
    fa.q = q;
    fa.cv = cv;
    fa.params = params;           % L (neuron x q), Ph, d - for fastfa_estep on new data
    fa.Lorth = Lorth;             % orthonormal loadings (= C)
    fa.TT = TT;                   % factors = TT * fastfa_estep(X, params).mean
    fa.explained = explained;     % % of shared variance per orthonormalised factor
    fa.pct_shared_variance = 100 * mean(shared);
    fa.traj = stack_traj(newD);   % trial x factor x time
    fa.mean_traj = squeeze(mean(fa.traj, 1));
    fa.loo_distance = loo_distance(fa.traj(:, 1:min(3, q), :));
    fa.C = C;
    pop.fa = fa;

    %% PCA
    n_pc = min(opts.n_pc, n);
    [newD, C, ~, explained] = reducedims_EH(D, 1, n_pc, []);
    pc.n_pc = n_pc;
    pc.coeff = C;                 % neuron x PC
    pc.mu = mean([D.data], 2);    % projections = coeff' * (X - mu)
    pc.explained = explained;     % % variance per PC (all PCs)
    pc.traj = stack_traj(newD);
    pc.mean_traj = squeeze(mean(pc.traj, 1));
    pc.loo_distance = loo_distance(pc.traj(:, 1:3, :));
    pop.pca = pc;

    pops(p) = pop; %#ok<AGROW>

    for method = ["PCA", "FA"]
        figs(end + 1) = plot_trajectories(pop, method, t, phase, trials, tag); %#ok<AGROW>
        names{end + 1} = sprintf('%s_%s_trajectories_%s', base, method, P.name); %#ok<AGROW>
        figs(end + 1) = plot_single_trials(pop, method, phase, tag); %#ok<AGROW>
        names{end + 1} = sprintf('%s_%s_single_trials_%s', base, method, P.name); %#ok<AGROW>
    end
end

results.subject = subject;
results.date = date;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.trials = trials;
results.t = t;
results.phase = phase;
results.phase_names = ["pre", string(trials.element_names), "post"];
results.pop = pops;
if ~opts.save; return; end

safe_save(fullfile(out_dir, [base '_PCA_FA_results.mat']), struct('results', results), cfg);
save_vector_figures(figs, names, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function opts = fill_defaults(opts)
def.populations = struct('name', {'CA1_good', 'CA1_all', 'V1_good'}, 'region', {"CA1", "CA1", "V1"}, ...
    'good_only', {true, false, true});
def.window = [];              % s around A onset; [] = [-0.5, sequence + 0.5] (-0.5 to 1.1 s for 150 ms)
def.bin = 0.02;               % s
def.sigma = 0.04;             % s, Gaussian smoothing
def.min_rate = 0.5;           % Hz, mean over trial windows
def.softnorm_const = 5;       % Hz (SOFTNORM_CONST in FA_continuous_test)
def.k_folds = 5;
def.max_q = 50;
def.n_pc = 10;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
end


function sm = smooth_counts(counts, opts)
sig = opts.sigma / opts.bin;
k = exp(-0.5 * ((-ceil(3 * sig):ceil(3 * sig)) / sig).^2);
k = reshape(k / sum(k), 1, 1, []);
sm = convn(counts, k, 'same');
end


function ph = phase_of(t, trials)
% 1 = before A, 2..5 = A..D (each element until the next one starts, i.e. including any
% grey between elements), 6 = after D
n = numel(trials.element_names);
ph = ones(size(t));
for e = 1:n
    if e < n; stop = trials.element_start(e + 1); else; stop = trials.element_end(e); end
    ph(t >= trials.element_start(e) & t < stop) = e + 1;
end
ph(t >= trials.element_end(n)) = n + 2;
end


function cv = fa_cross_validate(D, opts, n)
% 5-fold CV over trials as in FA_continuous_test: summed held-out log-likelihood
% (fastfa_estep) and leave-neuron-out squared error (cosmoother_fa) per q.
cv.candidate_dims = 5:min(n - 5, opts.max_q);
assert(~isempty(cv.candidate_dims), 'Too few neurons (%d) for FA cross-validation', n);
rng(opts.seed, 'twister');
n_tr = numel(D);
perm = randperm(n_tr);
fold = zeros(1, n_tr);
fold(perm) = mod(0:n_tr - 1, opts.k_folds) + 1;
cv.fold = fold;
cv.like = nan(numel(cv.candidate_dims), opts.k_folds);
cv.mse = nan(numel(cv.candidate_dims), opts.k_folds);
t0 = tic;
for i = 1:numel(cv.candidate_dims)
    for k = 1:opts.k_folds
        train = [D(fold ~= k).data];
        test = [D(fold == k).data];
        prm = fastfa(train, cv.candidate_dims(i), 'typ', 'fa');
        assert(~isempty(prm), 'FA CV failed at q = %d, fold %d', cv.candidate_dims(i), k);
        [~, cv.like(i, k)] = fastfa_estep(test, prm);
        cv.mse(i, k) = sum(sum((cosmoother_fa(test, prm) - test).^2));
    end
end
fprintf('FA cross-validation: %d dims x %d folds in %.0f s\n', numel(cv.candidate_dims), opts.k_folds, toc(t0));
cv.like_sum = sum(cv.like, 2);
cv.mse_sum = sum(cv.mse, 2);
[~, i_opt] = max(cv.like_sum);
cv.q_opt = cv.candidate_dims(i_opt);
end


function X = stack_traj(newD)
X = permute(cat(3, newD.data), [3 1 2]);   % trial x dim x time
end


function d = loo_distance(X)
% Euclidean distance of each trial to the mean of all other trials: trial x time
n = size(X, 1);
total = sum(X, 1);
d = zeros(n, size(X, 3));
for i = 1:n
    loo = (total - X(i, :, :)) / (n - 1);
    d(i, :) = squeeze(sqrt(sum((X(i, :, :) - loo).^2, 2)));
end
end


%% ------------------------------------------------------------------------
%% Plotting
function c = phase_colors()
% pre (grey), A, B, C, D (categorical slots blue, orange, aqua, violet), post (light grey)
c = [0.55 0.55 0.55; 0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655; 0.78 0.78 0.78];
end


function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end


function M = method_struct(pop, method)
if method == "PCA"
    M = pop.pca;
    M.label = 'PC';
    M.var_text = arrayfun(@(k) sprintf('%.0f%%', M.explained(k)), 1:3, 'UniformOutput', false);
else
    M = pop.fa;
    M.label = 'Factor ';
    M.var_text = arrayfun(@(k) sprintf('%.0f%% shared', M.explained(k)), 1:3, 'UniformOutput', false);
end
end


function draw_phases(ax, X, dims, phase, alpha, width)
% X: trial x dim x time. One line object per phase (trials joined with NaN gaps);
% each segment includes the next point so phases join up.
c = phase_colors();
for ph = 1:max(phase)
    idx = find(phase == ph);
    if isempty(idx); continue; end
    idx = [idx idx(end) + 1]; %#ok<AGROW>
    idx = idx(idx <= size(X, 3));
    seg = X(:, dims, idx);                                     % trial x d x t
    seg = cat(3, seg, nan(size(seg, 1), numel(dims), 1));
    xyz = reshape(permute(seg, [3 1 2]), [], numel(dims));     % (t+1)*trial x d
    col = [c(ph, :) alpha];
    if numel(dims) == 3
        plot3(ax, xyz(:, 1), xyz(:, 2), xyz(:, 3), '-', 'Color', col, 'LineWidth', width);
    else
        plot(ax, xyz(:, 1), xyz(:, 2), '-', 'Color', col, 'LineWidth', width);
    end
end
end


function shade_elements(ax, trials, yl)
c = phase_colors();
for e = 1:numel(trials.element_names)
    x0 = 1000 * trials.element_start(e);
    x1 = x0 + 1000 * trials.element_duration;
    patch(ax, [x0 x1 x1 x0], yl([1 1 2 2]), c(e + 1, :), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    text(ax, (x0 + x1) / 2, yl(2), trials.element_names{e}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
end


function f = plot_trajectories(pop, method, t, phase, trials, tag)
M =method_struct(pop, method);
X = M.traj;
mu = permute(M.mean_traj, [3 1 2]);   % 1 x dim x time
n_tr = size(X, 1);
f = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(f, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - %s, %s (n = %d units, %d sequences): thin = single trials, bold = mean', ...
    tag, method, strrep(pop.name, '_', ' '), numel(pop.cluster_id), n_tr), 'FontSize', 11);
lab = @(k) sprintf('%s%d (%s)', M.label, k, M.var_text{k});

% 3D: all trials + mean
ax = nexttile(tl, 1, [2 2]); hold(ax, 'on');
draw_phases(ax, X, 1:3, phase, 0.12, 0.5);
draw_phases(ax, mu, 1:3, phase, 1, 3.5);
[~, i0] = min(abs(t));
plot3(ax, mu(1, 1, i0), mu(1, 2, i0), mu(1, 3, i0), 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel(ax, lab(1)); ylabel(ax, lab(2)); zlabel(ax, lab(3));
view(ax, 45, 25); grid(ax, 'on'); axis(ax, 'tight');
c = phase_colors();
names = ["before A", string(trials.element_names), "after D"];
h = gobjects(1, 6);
for k = 1:6
    h(k) = plot3(ax, nan, nan, nan, '-', 'Color', c(k, :), 'LineWidth', 3);
end
legend(ax, h, names, 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 8);
title(ax, '3D trajectory (triangle = A onset)', 'FontWeight', 'normal');
style(ax);

% dims vs time
tiles = [3 4 7];
for k = 1:3
    ax = nexttile(tl, tiles(k)); hold(ax, 'on');
    Y = squeeze(X(:, k, :));
    lo = prctile(Y(:), 1); hi = prctile(Y(:), 99);
    yl = [lo hi] + [-0.1 0.25] * (hi - lo);
    shade_elements(ax, trials, yl);
    plot(ax, 1000 * t, Y', '-', 'Color', [0.2 0.2 0.2 0.06], 'LineWidth', 0.5);
    plot(ax, 1000 * t, squeeze(mu(1, k, :)), 'k-', 'LineWidth', 2.5);
    ylim(ax, yl); xlim(ax, 1000 * t([1 end]));
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, lab(k));
    style(ax);
end

% trial-to-trial variability
ax = nexttile(tl, 8); hold(ax, 'on');
d = M.loo_distance;
m = mean(d, 1);
s = std(d, 0, 1) / sqrt(size(d, 1));
yl = [0 1.25 * max(m + s)];
shade_elements(ax, trials, yl);
fill(ax, 1000 * [t fliplr(t)], [m - s fliplr(m + s)], [0.2 0.2 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ax, 1000 * t, m, 'k-', 'LineWidth', 2);
ylim(ax, yl); xlim(ax, 1000 * t([1 end]));
xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'distance to LOO mean (dims 1-3)');
title(ax, 'single-trial deviation from mean (\pm SEM)', 'FontWeight', 'normal');
style(ax);

% 2D: all trials, then the mean alone (zoomed; dots every 100 ms)
ax = nexttile(tl, 9); hold(ax, 'on');
draw_phases(ax, X, [1 2], phase, 0.12, 0.5);
draw_phases(ax, mu, [1 2], phase, 1, 3);
plot(ax, mu(1, 1, i0), mu(1, 2, i0), 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
xlabel(ax, lab(1)); ylabel(ax, lab(2)); axis(ax, 'tight');
title(ax, 'all trials + mean', 'FontWeight', 'normal');
style(ax);
tick = find(diff(floor(t / 0.1)) ~= 0) + 1;   % first sample after each 100 ms boundary
ax = nexttile(tl, 10); hold(ax, 'on');
draw_phases(ax, mu, 1:3, phase, 1, 3);
plot3(ax, squeeze(mu(1, 1, tick)), squeeze(mu(1, 2, tick)), squeeze(mu(1, 3, tick)), 'k.', 'MarkerSize', 10);
plot3(ax, mu(1, 1, i0), mu(1, 2, i0), mu(1, 3, i0), 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
xlabel(ax, lab(1)); ylabel(ax, lab(2)); zlabel(ax, lab(3));
view(ax, 45, 25); grid(ax, 'on'); axis(ax, 'tight');
title(ax, 'mean trajectory only (dots = 100 ms)', 'FontWeight', 'normal');
style(ax);
ax = nexttile(tl, 11); hold(ax, 'on');
draw_phases(ax, mu, [1 2], phase, 1, 3);
plot(ax, squeeze(mu(1, 1, tick)), squeeze(mu(1, 2, tick)), 'k.', 'MarkerSize', 10);
plot(ax, mu(1, 1, i0), mu(1, 2, i0), 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
xlabel(ax, lab(1)); ylabel(ax, lab(2)); axis(ax, 'tight');
title(ax, 'mean trajectory only', 'FontWeight', 'normal');
style(ax);

% variance / dimensionality
ax = nexttile(tl, 12); hold(ax, 'on');
if method == "PCA"
    ev = M.explained(1:min(20, end));
    bar(ax, ev, 'FaceColor', [0.165 0.471 0.839], 'EdgeColor', 'none');
    plot(ax, cumsum(ev), 'k-o', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    xlabel(ax, 'PC'); ylabel(ax, '% variance (bars), cumulative (line)');
    title(ax, sprintf('PCA: top 3 PCs = %.0f%% of variance', sum(M.explained(1:3))), 'FontWeight', 'normal');
else
    ll = M.cv.like_sum;
    plot(ax, M.cv.candidate_dims, ll - max(ll), 'k-o', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    xline(ax, M.q, '-', sprintf('q = %d', M.q), 'Color', [0.922 0.408 0.204], 'LineWidth', 1.5, 'FontSize', 8);
    xlabel(ax, 'number of factors q'); ylabel(ax, 'held-out log-likelihood (re max)');
    title(ax, sprintf('FA %d-fold CV over trials; shared variance %.0f%%', size(M.cv.like, 2), ...
        M.pct_shared_variance), 'FontWeight', 'normal');
end
style(ax);
end


function f = plot_single_trials(pop, method, phase, tag)
% Every sequence as its own small panel (dims 1 vs 2), mean trajectory in grey.
M = method_struct(pop, method);
X = M.traj;
mu = permute(M.mean_traj, [3 1 2]);
n_tr = size(X, 1);
n_col = 20;
n_row = ceil(n_tr / n_col);
f = figure('Color', 'w', 'Position', [20 20 2200 1250]);
tl = tiledlayout(f, n_row, n_col, 'TileSpacing', 'none', 'Padding', 'compact');
title(tl, sprintf(['%s - %s, %s: every sequence, %s1 (x) vs %s2 (y), same axes; grey = mean trajectory; ' ...
    'colour = before A / A / B / C / D / after D'], tag, method, strrep(pop.name, '_', ' '), M.label, M.label), ...
    'FontSize', 11);
xl = [prctile(X(:, 1, :), 0.5, 'all') prctile(X(:, 1, :), 99.5, 'all')];
yl = [prctile(X(:, 2, :), 0.5, 'all') prctile(X(:, 2, :), 99.5, 'all')];
for i = 1:n_tr
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, squeeze(mu(1, 1, :)), squeeze(mu(1, 2, :)), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
    draw_phases(ax, X(i, :, :), [1 2], phase, 1, 1);
    xlim(ax, xl); ylim(ax, yl);
    set(ax, 'XTick', [], 'YTick', [], 'Box', 'on', 'XColor', [0.85 0.85 0.85], 'YColor', [0.85 0.85 0.85]);
    text(ax, xl(1), yl(2), sprintf(' %d', i), 'VerticalAlignment', 'top', 'FontSize', 6, 'Color', [0.4 0.4 0.4]);
end
end
