function results = run_abcd_seqlik_contributions(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_SEQLIK_CONTRIBUTIONS  Which neurons carry the sequence-likelihood score?
%
%   results = run_abcd_seqlik_contributions('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Uses the folds, units and grey windows of run_abcd_sequence_likelihood (run it first).
% Every held-out sequence (and grey window) is re-scored with the matrix of its own CV
% fold, and Lseq is split into per-neuron contributions (seq_likelihood, do_contrib):
% each spike-pair term goes half to each of its two neurons, minus the neuron's mean
% contribution in the cell-ID shuffles, divided by the shuffle SD. A neuron's
% contribution is therefore in z units, and the contributions of a window add up to
% its z. Positive = this neuron's spikes sit in the trained order more than chance.
%
% Per population: mean contribution per neuron over held-out sequences (signrank vs 0,
% FDR), how concentrated the total is (share carried by the top neurons; effective
% number of neurons per sequence, N_eff = (sum |c|)^2 / sum c^2), and how contribution
% relates to firing rate and to the neuron's order strength (mean |W| of its row).
%
% Outputs (decoding_* names) in <analysis>/<subject>/<date>/<stimulus>/Decoding.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.populations = {'CA1_all', 'CA1_good', 'V1_good'};
def.n_shuffle = 500;
def.fdr = 0.05;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);

paths  = get_session_paths(subject, date, stimulus_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
S = load(fullfile(out_dir, sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
SL = S.results;
trials = SL.trials;
win = SL.opts.window;
fold = SL.fold;
n_seq = numel(trials.seq_onset);
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));

for p = 1:numel(opts.populations)
    name = opts.populations{p};
    P = SL.pop(strcmp({SL.pop.name}, name));
    ids = P.cluster_id;
    n = numel(ids);
    stim_w = collect_spike_windows(units.spike_times, units.spike_id, ids, trials.seq_onset + win(1), trials.seq_onset + win(2));
    grey_w = collect_spike_windows(units.spike_times, units.spike_id, ids, SL.grey_start, SL.grey_start + diff(win));
    Cs = nan(n, n_seq); Cg = nan(n, n_seq); zs = nan(1, n_seq); zg = nan(1, n_seq);
    t0 = tic;
    for f = 1:max(fold)
        M = seq_order_matrix(stim_w(fold ~= f), n);
        W = M.W;
        test = find(fold == f);
        ns = opts.n_shuffle;
        cs = cell(numel(test), 1); cg = cell(numel(test), 1); a = zeros(numel(test), 1); b = a;
        parfor k = 1:numel(test)
            sw = stim_w{test(k)}; gw = grey_w{test(k)}; %#ok<PFBNS>
            r = seq_likelihood(sw(:, 1), sw(:, 2), W, ns, test(k), true);
            cs{k} = r.contrib_z; a(k) = r.z;
            r = seq_likelihood(gw(:, 1), gw(:, 2), W, ns, n_seq + test(k), true);
            cg{k} = r.contrib_z; b(k) = r.z;
        end
        Cs(:, test) = [cs{:}]; Cg(:, test) = [cg{:}]; zs(test) = a; zg(test) = b;
    end
    fprintf('%s: %d neurons, decomposed %d + %d windows in %.0f s (check: max |sum(contrib) - z| = %.1e)\n', ...
        name, n, n_seq, n_seq, toc(t0), max(abs(sum(Cs, 1) - zs)));

    M_all = seq_order_matrix(stim_w, n);
    R.name = name;
    R.cluster_id = ids;
    [~, row] = ismember(ids, units.cluster_id);
    R.depth = units.depth(row);
    R.shank = units.shank(row);
    R.is_good = units.is_good(row);
    R.cell_type = units.cell_type(row);
    R.rate_hz = accum_rate(stim_w, n, diff(win));
    R.order_index = (sum(M_all.p, 2) - 0.5) / (n - 1);
    W0 = M_all.W; W0(1:n + 1:end) = NaN;
    R.order_strength = mean(abs(W0), 2, 'omitnan');           % how reliable this neuron's order is
    R.contrib_stim = Cs;                                       % neuron x sequence (z units)
    R.contrib_grey = Cg;
    R.z_stim = zs; R.z_grey = zg;
    R.mean_stim = mean(Cs, 2, 'omitnan');
    R.sem_stim = std(Cs, 0, 2, 'omitnan') / sqrt(n_seq);
    R.mean_grey = mean(Cg, 2, 'omitnan');
    R.frac_positive = mean(Cs > 0, 2);
    R.p = arrayfun(@(i) signrank(Cs(i, :)), (1:n)');
    R.sig = fdr_bh(R.p, opts.fdr);
    pos = max(R.mean_stim, 0);
    srt = sort(pos, 'descend');
    R.cum_share = cumsum(srt) / sum(srt);                       % share of the positive total by top k
    R.n_for_50 = find(R.cum_share >= 0.5, 1);
    R.n_for_80 = find(R.cum_share >= 0.8, 1);
    R.n_eff = (sum(abs(Cs), 1).^2) ./ sum(Cs.^2, 1);           % per sequence
    srt_tr = sort(Cs, 1, 'descend');
    R.top3_share = sum(srt_tr(1:min(3, n), :), 1) ./ sum(Cs, 1);
    R.top3_share(zs <= 0) = NaN;                                % only meaningful for forward-scoring sequences
    fprintf(['   significant neurons (FDR %.2f): %d positive, %d negative | top neurons for 50%% / 80%% of the ' ...
        'positive total: %d / %d of %d | median N_eff per sequence %.1f | median top-3 share %.0f%%\n'], ...
        opts.fdr, sum(R.sig & R.mean_stim > 0), sum(R.sig & R.mean_stim < 0), R.n_for_50, R.n_for_80, n, ...
        median(R.n_eff, 'omitnan'), 100 * median(R.top3_share, 'omitnan'));
    pops(p) = R; %#ok<AGROW>
end

results.subject = subject;
results.date = date;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.pop = pops;
figs = [plot_contrib(pops, opts, tag), plot_heatmaps(pops, tag)];
names = {sprintf('decoding_%s_seqlik_contributions', stimulus_name), sprintf('decoding_%s_seqlik_contrib_heatmap', stimulus_name)};
if ~opts.save; return; end
safe_save(fullfile(out_dir, sprintf('decoding_%s_seqlik_contrib_results.mat', stimulus_name)), struct('results', results), cfg);
save_vector_figures(figs, names, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function r = accum_rate(w, n, dur)
c = zeros(n, 1);
for k = 1:numel(w)
    if ~isempty(w{k}); c = c + accumarray(w{k}(:, 2), 1, [n 1]); end
end
r = c / (numel(w) * dur);
end

function sig = fdr_bh(p, q)
[ps, o] = sort(p(:));
m = numel(ps);
k = find(ps <= (1:m)' / m * q, 1, 'last');
sig = false(size(p));
if ~isempty(k); sig(o(1:k)) = true; end
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function cmap = div_cmap()
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 0.94 0.94 0.94; 0.922 0.408 0.204], linspace(0, 1, 256));
end


function f = plot_contrib(pops, opts, tag)
np = numel(pops);
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.6 0.6 0.6];
f = figure('Color', 'w', 'Position', [30 30 1700 330 * np]);
tl = tiledlayout(f, np, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - which neurons carry L_{seq}? Per-neuron contribution (z units; a sequence''s contributions ' ...
    'add up to its z), held-out sequences, %d shuffles'], tag, opts.n_shuffle), 'FontSize', 11);
for p = 1:np
    R = pops(p);
    n = numel(R.cluster_id);
    pname = sprintf('%s (n = %d)', strrep(R.name, '_', ' '), n);

    % 1. per-neuron mean contribution, sorted
    ax = nexttile(tl); hold(ax, 'on');
    [~, o] = sort(R.mean_stim, 'descend');
    x = 1:n;
    b = bar(ax, x, R.mean_stim(o), 1, 'FaceColor', 'flat', 'EdgeColor', 'none');
    col = repmat(grey, n, 1);
    col(R.sig(o) & R.mean_stim(o) > 0, :) = repmat(orange, sum(R.sig(o) & R.mean_stim(o) > 0), 1);
    col(R.sig(o) & R.mean_stim(o) < 0, :) = repmat(blue, sum(R.sig(o) & R.mean_stim(o) < 0), 1);
    b.CData = col;
    errorbar(ax, x, R.mean_stim(o), R.sem_stim(o), 'k', 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.5);
    plot(ax, x, R.mean_grey(o), '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 6);
    yline(ax, 0, 'k-');
    xlim(ax, [0 n + 1]);
    xlabel(ax, 'neuron (sorted by mean contribution)'); ylabel(ax, 'mean contribution (z units)');
    title(ax, sprintf('%s: %d neurons significantly + (orange), %d - (blue); dots = grey windows', pname, ...
        sum(R.sig & R.mean_stim > 0), sum(R.sig & R.mean_stim < 0)), 'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 2. concentration
    ax = nexttile(tl); hold(ax, 'on');
    k = (1:n)';
    plot(ax, [0; k], [0; R.cum_share], '-', 'Color', orange, 'LineWidth', 2);
    plot(ax, [0 n], [0 1], 'k:');
    plot(ax, R.n_for_50, 0.5, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    plot(ax, R.n_for_80, 0.8, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
    xlim(ax, [0 n]); ylim(ax, [0 1.02]);
    xlabel(ax, 'top-k neurons'); ylabel(ax, 'share of the positive total');
    title(ax, sprintf('50%% from %d neurons, 80%% from %d (of %d); dotted = equal shares', R.n_for_50, R.n_for_80, n), ...
        'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 3. per-sequence concentration
    ax = nexttile(tl); hold(ax, 'on');
    histogram(ax, R.n_eff, 25, 'FaceColor', blue, 'EdgeColor', 'none');
    xline(ax, median(R.n_eff, 'omitnan'), 'k-', sprintf('median %.1f', median(R.n_eff, 'omitnan')), 'FontSize', 7);
    xlabel(ax, 'N_{eff} per held-out sequence = (\Sigma|c|)^2 / \Sigmac^2'); ylabel(ax, 'sequences');
    title(ax, sprintf('effective number of contributing neurons; top-3 share median %.0f%% (z > 0)', ...
        100 * median(R.top3_share, 'omitnan')), 'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 4. contribution vs firing rate, coloured by order index
    ax = nexttile(tl); hold(ax, 'on');
    scatter(ax, R.rate_hz, R.mean_stim, 18 + 30 * R.sig, R.order_index, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
    colormap(ax, div_cmap()); cb = colorbar(ax); cb.Label.String = 'order index (early \rightarrow late)';
    set(ax, 'XScale', 'log');
    yline(ax, 0, 'k:');
    rho = corr(R.rate_hz, R.mean_stim, 'type', 'Spearman');
    rho2 = corr(R.order_strength, R.mean_stim, 'type', 'Spearman');
    xlabel(ax, 'firing rate in sequence windows (Hz)'); ylabel(ax, 'mean contribution (z units)');
    title(ax, sprintf('vs rate \\rho = %.2f; vs order strength \\rho = %.2f (big dots = significant)', rho, rho2), ...
        'FontWeight', 'normal', 'FontSize', 8);
    style(ax);
end
end


function f = plot_heatmaps(pops, tag)
np = numel(pops);
f = figure('Color', 'w', 'Position', [30 30 520 * np 620]);
tl = tiledlayout(f, 1, np, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - per-neuron contribution in every held-out sequence (rows: neurons sorted by mean; ' ...
    'columns: sequences sorted by z, best on the left)'], tag), 'FontSize', 11);
for p = 1:np
    R = pops(p);
    [~, o] = sort(R.mean_stim, 'descend');
    [~, s] = sort(R.z_stim, 'descend');
    ax = nexttile(tl);
    C = R.contrib_stim(o, s);
    lim = max(0.3, prctile(abs(C(:)), 99));
    imagesc(ax, C, [-lim lim]); colormap(ax, div_cmap());
    cb = colorbar(ax); cb.Label.String = 'contribution (z units)';
    xlabel(ax, 'held-out sequence (sorted by z)'); ylabel(ax, 'neuron (sorted by mean contribution)');
    title(ax, sprintf('%s (n = %d)', strrep(R.name, '_', ' '), numel(R.cluster_id)), 'FontWeight', 'normal');
    style(ax);
end
end
