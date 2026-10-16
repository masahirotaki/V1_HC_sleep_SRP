function results = run_abcd_sequence_likelihood(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_SEQUENCE_LIKELIHOOD  Is there a consistent firing ORDER in CA1 / V1 during
% the GAVNIK ABCD sequence? Sequence-likelihood method of Huh, Yun, ..., Jung (2026)
% Nat Commun 17:8304 (seq_order_matrix / seq_likelihood), tested on awake data with
% cross-validation before it is applied to ripple events.
%
%   results = run_abcd_sequence_likelihood('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% 1. Active-phase windows: 0-0.8 s after A onset of each sequence (A-D = 0-0.6 s,
%    + CA1 response latency). All spikes of the population in the window.
% 2. 10-fold CV over sequences: the ordered firing probability matrix is built from
%    the training sequences, then every held-out sequence is scored (Lseq).
% 3. Two nulls:
%    - cell-identity shuffle of each held-out window (1000 shuffles; z-score and
%      percentile of Lseq), as in the paper;
%    - grey windows: one 0.8 s window per sequence at a random position in the grey
%      screen (>= 0.3 s after D offset, ending before the next A), scored with the same
%      fold's matrix and its own shuffle null.
% A sequence is significant forward if its Lseq exceeds 95% of its shuffles (one-tailed;
% chance 5%). Performance: z-score distributions, ROC AUC (stimulus vs grey), fraction
% significant, and fraction of sequences above the grey windows' 95th percentile z.
% NB: in awake 0.8 s windows the cell-ID shuffle null is too liberal (grey windows are
% "significant" well above 5% in BOTH tails, because spikes cluster in time - bursts,
% population events, slow rate changes - which the shuffle destroys). The grey windows
% are therefore the calibrated null; the shuffle z is used as a normalised score.
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
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));
if isempty(opts.window)   % A onset to D offset + 0.2 s (CA1 latency), nominal element duration
    opts.window = [0, round(trials.seq_duration, 2) + 0.2];
end
dur = diff(opts.window);

%% windows
rng(opts.seed);
stim_start = trials.seq_onset + opts.window(1);
% grey pool: [D offset + margin, next A onset) of every gap between sequences
g0 = trials.seq_offset(1:end - 1) + opts.grey_margin;
g1 = trials.seq_onset(2:end);
ok = g1 - g0 >= dur;
g0 = g0(ok); g1 = g1(ok);
span = g1 - g0 - dur;                                   % room for the window start
switch opts.grey_sampling
    case 'random'       % original: start uniform over all eligible grey (windows may overlap)
        pick = randsample(numel(span), n_seq, true, span / sum(span));
        grey_start = g0(pick) + rand(n_seq, 1) .* span(pick);
        grey_in_pause = (g1(pick) - g0(pick)) > opts.pause_min;
    case 'nonoverlap'   % one random window per short gap + non-overlapping tiles in pauses
        cand = []; is_pause = [];
        for k = 1:numel(g0)
            if g1(k) - g0(k) > opts.pause_min
                s = (g0(k):dur:g1(k) - dur)';
                cand = [cand; s]; is_pause = [is_pause; true(numel(s), 1)]; %#ok<AGROW>
            else
                cand = [cand; g0(k) + rand * span(k)]; is_pause = [is_pause; false]; %#ok<AGROW>
            end
        end
        pick = randsample(numel(cand), min(n_seq, numel(cand)));
        if numel(pick) < n_seq; pick = [pick; randsample(numel(cand), n_seq - numel(pick))]; end
        grey_start = cand(pick);
        grey_in_pause = is_pause(pick);
    otherwise
        error('Unknown grey_sampling %s', opts.grey_sampling);
end
fold = mod(randperm(n_seq) - 1, opts.n_folds) + 1;
fprintf('%d stimulus windows, %d grey windows (%.1f s of grey available)\n', n_seq, n_seq, sum(span + dur));

base = sprintf('decoding_%s', stimulus_name);
for p = 1:numel(opts.populations)
    P = opts.populations(p);
    sel = units.region == P.region & (units.is_good | ~P.good_only);
    ids = units.cluster_id(sel);
    [stim_w, n_sp] = collect_windows(units, ids, stim_start, dur);
    keep = n_sp >= opts.min_spikes;
    ids = ids(keep);
    stim_w = collect_windows(units, ids, stim_start, dur);
    grey_w = collect_windows(units, ids, grey_start, dur);
    n = numel(ids);
    fprintf('\n=== %s: %d units (>= %d spikes in the stimulus windows)\n', P.name, n, opts.min_spikes);

    S = struct('L', [], 'L_norm', [], 'z', [], 'pct', [], 'n_spikes', [], 'n_neurons', []);
    stim = repmat(S, n_seq, 1);
    grey = repmat(S, n_seq, 1);
    t0 = tic;
    for f = 1:opts.n_folds
        M = seq_order_matrix(stim_w(fold ~= f), n);
        test = find(fold == f);
        W = M.W;
        ns = opts.n_shuffle;
        rs = cell(numel(test), 1);
        rg = cell(numel(test), 1);
        parfor k = 1:numel(test)
            sw = stim_w{test(k)}; gw = grey_w{test(k)}; %#ok<PFBNS>
            rs{k} = slim(seq_likelihood(sw(:, 1), sw(:, 2), W, ns, test(k)));
            rg{k} = slim(seq_likelihood(gw(:, 1), gw(:, 2), W, ns, n_seq + test(k)));
        end
        stim(test) = [rs{:}];
        grey(test) = [rg{:}];
    end
    fprintf('scored %d + %d windows in %.0f s\n', n_seq, n_seq, toc(t0));

    pop = struct('name', P.name, 'region', P.region, 'good_only', P.good_only, 'cluster_id', ids);
    pop.stim = stim;
    pop.grey = grey;
    pop.order_all = seq_order_matrix(stim_w, n);           % on all sequences, for display / ripples
    pop.order_index = (sum(pop.order_all.p, 2) - 0.5) / (n - 1);  % mean P(i fires after j): late = high
    pop.psth = window_psth(stim_w, n, dur, opts.psth_bin); % neuron x time (Hz), mean over sequences
    pop.example_stim = stim_w;  %#ok<*STRNU>
    pop.example_grey = grey_w;
    pop.summary = summarise([stim.z], [grey.z], [stim.pct], [grey.pct], opts.alpha);
    s = pop.summary;
    fprintf(['z: stimulus median %.2f vs grey %.2f (ranksum p = %.2g), AUC %.3f | significant forward: ' ...
        'stimulus %.1f%% (binomial p = %.2g), grey %.1f%% (p = %.2g)\n'], s.median_z_stim, s.median_z_grey, ...
        s.p_ranksum, s.auc, 100 * s.frac_fwd_stim, s.p_binom_stim, 100 * s.frac_fwd_grey, s.p_binom_grey);
    pops(p) = pop; %#ok<AGROW>
end

results.subject = subject;
results.date = date;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.trials = trials;
results.fold = fold;
results.grey_start = grey_start;
results.grey_in_pause = grey_in_pause;
results.pop = pops;
results.summary = summary_table(pops);
disp(results.summary);

figs = [plot_performance(pops, opts, tag), plot_order(pops, opts, trials, tag)];
names = {[base '_seqlik_performance'], [base '_seqlik_order']};
if ~opts.save; return; end
for p = 1:numel(results.pop)   % spike windows are only needed for the figures
    results.pop(p).example_stim = [];
    results.pop(p).example_grey = [];
end
safe_save(fullfile(out_dir, [base '_seqlik_results.mat']), struct('results', results), cfg);
csv_file = fullfile(out_dir, [base '_seqlik_summary.csv']);
check_output_path(csv_file, cfg);
check_overwrite(csv_file, cfg);
writetable(results.summary, csv_file);
save_vector_figures(figs, names, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function opts = fill_defaults(opts)
def.populations = struct('name', {'CA1_good', 'CA1_all', 'V1_good'}, 'region', {"CA1", "CA1", "V1"}, ...
    'good_only', {true, false, true});
def.window = [];            % [] = [0, sequence + 0.2] s after A onset (0-0.8 s for 150 ms elements)
def.grey_margin = 0.3;      % s after D offset before grey windows may start (skip offset responses)
def.grey_sampling = 'random';   % 'random' (original, may overlap) or 'nonoverlap'
def.pause_min = 5;          % s; gaps longer than this are block pauses
def.n_folds = 10;
def.n_shuffle = 1000;
def.min_spikes = 20;        % per unit, summed over all stimulus windows
def.alpha = 0.05;           % one-tailed, forward
def.psth_bin = 0.01;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
end


function [w, n_sp] = collect_windows(units, ids, starts, dur)
% w{k} = [time re window start, neuron index] of all spikes of ids in window k
in = ismember(units.spike_id, ids);
st = units.spike_times(in);
[~, col] = ismember(units.spike_id(in), ids);
[st, o] = sort(st);
col = col(o);
w = cell(numel(starts), 1);
n_sp = zeros(numel(ids), 1);
for k = 1:numel(starts)
    a = find(st >= starts(k), 1);
    b = find(st < starts(k) + dur, 1, 'last');
    if isempty(a) || isempty(b) || b < a
        w{k} = zeros(0, 2);
        continue
    end
    w{k} = [st(a:b) - starts(k), col(a:b)];
    n_sp = n_sp + accumarray(col(a:b), 1, [numel(ids) 1]);
end
end


function r = slim(R)
r = struct('L', R.L, 'L_norm', R.L_norm, 'z', R.z, 'pct', R.pct, 'n_spikes', R.n_spikes, 'n_neurons', R.n_neurons);
end


function P = window_psth(w, n, dur, bin)
edges = 0:bin:dur + 1e-9;
P = zeros(n, numel(edges) - 1);
for k = 1:numel(w)
    if isempty(w{k}); continue; end
    b = discretize(w{k}(:, 1), edges);
    ok = ~isnan(b);
    P = P + accumarray([w{k}(ok, 2) b(ok)], 1, size(P));
end
P = P / (numel(w) * bin);
end


function s = summarise(zs, zg, ps, pg, alpha)
s.median_z_stim = median(zs, 'omitnan');
s.median_z_grey = median(zg, 'omitnan');
s.mean_z_stim = mean(zs, 'omitnan');
s.mean_z_grey = mean(zg, 'omitnan');
s.p_ranksum = ranksum(zs(~isnan(zs)), zg(~isnan(zg)));
s.p_signrank_stim = signrank(zs(~isnan(zs)));
s.p_signrank_grey = signrank(zg(~isnan(zg)));
[~, ~, ~, s.auc] = perfcurve([ones(1, sum(~isnan(zs))) zeros(1, sum(~isnan(zg)))], [zs(~isnan(zs)) zg(~isnan(zg))], 1);
s.frac_fwd_stim = mean(ps > 1 - alpha);
s.frac_fwd_grey = mean(pg > 1 - alpha);
s.frac_rev_stim = mean(ps < alpha);
s.frac_rev_grey = mean(pg < alpha);
s.p_binom_stim = 1 - binocdf(sum(ps > 1 - alpha) - 1, numel(ps), alpha);
s.p_binom_grey = 1 - binocdf(sum(pg > 1 - alpha) - 1, numel(pg), alpha);
% calibrated against the grey windows instead of the cell-ID shuffle (which is too
% liberal for long awake windows: grey windows exceed it in both tails)
s.grey_z95 = prctile(zg, 100 * (1 - alpha));
s.frac_above_grey95 = mean(zs > s.grey_z95);
s.p_binom_above_grey95 = 1 - binocdf(sum(zs > s.grey_z95) - 1, sum(~isnan(zs)), alpha);
end


function T = summary_table(pops)
rows = {};
for p = 1:numel(pops)
    s = pops(p).summary;
    rows(end + 1, :) = {pops(p).name, numel(pops(p).cluster_id), s.median_z_stim, s.median_z_grey, s.p_ranksum, ...
        s.auc, s.frac_fwd_stim, s.p_binom_stim, s.frac_fwd_grey, s.p_binom_grey, s.frac_rev_stim, s.frac_rev_grey, ...
        s.grey_z95, s.frac_above_grey95, s.p_binom_above_grey95}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'population', 'n_units', 'median_z_stimulus', 'median_z_grey', ...
    'p_ranksum', 'AUC_stimulus_vs_grey', 'frac_forward_stimulus', 'p_binom_stimulus', 'frac_forward_grey', ...
    'p_binom_grey', 'frac_reverse_stimulus', 'frac_reverse_grey', 'grey_z_95pct', 'frac_stimulus_above_grey95', ...
    'p_binom_above_grey95'});
end


%% ------------------------------------------------------------------------
%% Plotting
function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function c = cols()
c.stim = [0.165 0.471 0.839];
c.grey = [0.55 0.55 0.55];
end

function cmap = div_cmap()
% blue (earlier) - white (0.5) - orange (later)
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 1 1 1; 0.922 0.408 0.204], linspace(0, 1, 256));
end

function cmap = seq_cmap()
cmap = interp1([0 0.6 1], [1 1 1; 0.165 0.471 0.839; 0.04 0.15 0.35], linspace(0, 1, 256));
end


function f = plot_performance(pops, opts, tag)
np = numel(pops);
c = cols();
f = figure('Color', 'w', 'Position', [40 40 1600 300 * np]);
tl = tiledlayout(f, np, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - sequence likelihood (Huh et al. 2026): order matrix from training sequences ' ...
    '(%g-%g s after A), held-out sequences vs random grey windows, %d-fold CV, %d cell-ID shuffles'], ...
    tag, opts.window, opts.n_folds, opts.n_shuffle), 'FontSize', 11);
for p = 1:np
    P = pops(p);
    zs = [P.stim.z]; zg = [P.grey.z];
    s = P.summary;
    pname = sprintf('%s (n = %d)', strrep(P.name, '_', ' '), numel(P.cluster_id));

    ax = nexttile(tl); hold(ax, 'on');
    edges = linspace(min([zs zg -4]), max([zs zg 4]), 40);
    histogram(ax, zg, edges, 'FaceColor', c.grey, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'grey windows');
    histogram(ax, zs, edges, 'FaceColor', c.stim, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'held-out sequences');
    xline(ax, 0, 'k:', 'HandleVisibility', 'off');
    xline(ax, norminv(1 - opts.alpha), '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, s.grey_z95, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', ...
        sprintf('grey 95th pct: %.0f%% of sequences above', 100 * s.frac_above_grey95));
    xlabel(ax, 'L_{seq} z-score vs cell-ID shuffle'); ylabel(ax, 'windows');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s: median z %.2f vs %.2f (p = %.2g)', pname, s.median_z_stim, s.median_z_grey, s.p_ranksum), ...
        'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    ok_s = ~isnan(zs); ok_g = ~isnan(zg);
    [fx, fy] = perfcurve([ones(1, sum(ok_s)) zeros(1, sum(ok_g))], [zs(ok_s) zg(ok_g)], 1);
    plot(ax, fx, fy, '-', 'Color', c.stim, 'LineWidth', 2);
    plot(ax, [0 1], [0 1], 'k:');
    axis(ax, 'square'); xlabel(ax, 'false positive rate (grey)'); ylabel(ax, 'true positive rate (sequences)');
    title(ax, sprintf('stimulus vs grey: AUC = %.3f', s.auc), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    vals = 100 * [s.frac_fwd_stim s.frac_fwd_grey; s.frac_rev_stim s.frac_rev_grey];
    b = bar(ax, vals, 'BarWidth', 0.8, 'EdgeColor', 'none');
    b(1).FaceColor = c.stim; b(2).FaceColor = c.grey;
    % 95% upper bound of the binomial count at chance
    lim = 100 * binoinv(0.95, numel(zs), opts.alpha) / numel(zs);
    yline(ax, 100 * opts.alpha, ':', 'chance', 'Color', 'k', 'FontSize', 7);
    yline(ax, lim, '--', 'binomial 95%', 'Color', [0.3 0.3 0.3], 'FontSize', 7);
    for k = 1:2
        for m = 1:2
            text(ax, b(m).XEndPoints(k), vals(k, m) + 1, sprintf('%.0f%%', vals(k, m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 7);
        end
    end
    set(ax, 'XTick', 1:2, 'XTickLabel', {'forward (p > 0.95)', 'reverse (p < 0.05)'});
    ylabel(ax, '% of windows significant'); ylim(ax, [0 max(30, max(vals(:)) + 8)]);
    legend(ax, {'held-out sequences', 'grey windows'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('significant: %.0f%% vs %.0f%% (binomial p = %.2g)', 100 * s.frac_fwd_stim, ...
        100 * s.frac_fwd_grey, s.p_binom_stim), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    x = 1:numel(zs);
    plot(ax, x, zg, '.', 'Color', c.grey, 'MarkerSize', 8, 'DisplayName', 'grey');
    plot(ax, x, zs, '.', 'Color', c.stim, 'MarkerSize', 8, 'DisplayName', 'sequence');
    plot(ax, x, movmean(zs, 15, 'omitnan'), '-', 'Color', c.stim, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ax, x, movmean(zg, 15, 'omitnan'), '-', 'Color', c.grey, 'LineWidth', 2, 'HandleVisibility', 'off');
    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    xlabel(ax, 'sequence # (session order)'); ylabel(ax, 'L_{seq} z');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, 'across the session (line = 15-sequence mean)', 'FontWeight', 'normal');
    style(ax);
end
end


function f = plot_order(pops, opts, trials, tag)
np = numel(pops);
f = figure('Color', 'w', 'Position', [40 40 1600 330 * np]);
tl = tiledlayout(f, np, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - ordered firing probability matrix (all sequences, %g-%g s) and example windows; ' ...
    'neurons sorted by order index (mean P(fires after others))'], tag, opts.window), 'FontSize', 11);
el = 1000 * unique([trials.element_start(2:end) trials.element_end]);
for p = 1:np
    P = pops(p);
    [~, o] = sort(P.order_index);
    n = numel(o);
    pname = sprintf('%s (n = %d)', strrep(P.name, '_', ' '), n);

    ax = nexttile(tl);
    v = P.order_all.p(o, o);
    lim = max(0.02, prctile(abs(v(:) - 0.5), 98));
    imagesc(ax, v, 0.5 + [-lim lim]); colormap(ax, div_cmap()); axis(ax, 'square'); set(ax, 'YDir', 'normal');
    cb = colorbar(ax); cb.Label.String = 'P(row fires after column)';
    xlabel(ax, 'neuron (sorted, early \rightarrow late)'); ylabel(ax, 'neuron (sorted)');
    title(ax, [pname ': order matrix'], 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl);
    R = P.psth(o, :);
    R = (R - min(R, [], 2)) ./ max(max(R, [], 2) - min(R, [], 2), eps);
    t = 1000 * ((1:size(R, 2)) - 0.5) * opts.psth_bin + 1000 * opts.window(1);
    imagesc(ax, t, 1:n, movmean(R, 3, 2), [0 1]); colormap(ax, seq_cmap()); set(ax, 'YDir', 'normal');
    xline(ax, [0 el], 'k:');
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'neuron (same order)');
    title(ax, 'mean response, same sort (min-max per neuron)', 'FontWeight', 'normal');
    style(ax);

    [~, best] = max([P.stim.z]);
    raster(nexttile(tl), P.example_stim{best}, o, opts, el, ...
        sprintf('held-out sequence %d: z = %.1f (best)', best, P.stim(best).z), [0.165 0.471 0.839]);
    zg = [P.grey.z];
    [~, med] = min(abs(zg - median(zg, 'omitnan')));
    raster(nexttile(tl), P.example_grey{med}, o, opts, [], ...
        sprintf('grey window (median z): z = %.1f', P.grey(med).z), [0.4 0.4 0.4]);
end
end


function raster(ax, w, o, opts, el, ttl, col)
hold(ax, 'on');
rank = zeros(numel(o), 1);
rank(o) = 1:numel(o);
if ~isempty(w)
    plot(ax, 1000 * (w(:, 1) + opts.window(1)), rank(w(:, 2)), '.', 'Color', col, 'MarkerSize', 5);
end
if ~isempty(el); xline(ax, [0 el], 'k:'); end
xlim(ax, 1000 * opts.window); ylim(ax, [0 numel(o) + 1]);
xlabel(ax, 'time in window (ms)'); ylabel(ax, 'neuron (sorted early \rightarrow late)');
title(ax, ttl, 'FontWeight', 'normal');
style(ax);
end
