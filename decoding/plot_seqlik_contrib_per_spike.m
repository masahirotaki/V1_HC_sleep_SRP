function plot_seqlik_contrib_per_spike(subject, date, stimulus_name, cfg)
% PLOT_SEQLIK_CONTRIB_PER_SPIKE  Per-spike view of the Lseq contributions from
% run_abcd_seqlik_contributions.
%
%   plot_seqlik_contrib_per_spike('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% A neuron's contribution to a held-out sequence (z units) = spikes it fired in that
% window x contribution per spike. Per spike = c(i, t) / n(i, t) for windows where the
% neuron fired; per neuron: mean over those windows +- SEM, signrank vs 0 (FDR 0.05).
% Panels per population: per-spike contribution ranked; the decomposition (mean spikes per
% window x per spike = mean total, iso-lines); per spike vs order index and vs order
% strength; and a neuron x sequence heatmap of per-spike contribution.
% Saved as decoding_<stimulus>_seqlik_contrib_per_spike(_heatmap) in the Decoding folder.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths  = get_session_paths(subject, date, stimulus_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
S = load(fullfile(out_dir, sprintf('decoding_%s_seqlik_contrib_results.mat', stimulus_name)));
C = S.results;
L = load(fullfile(out_dir, sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
tr = L.results.trials;
win = L.results.opts.window;
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));

for p = 1:numel(C.pop)
    R = C.pop(p);
    n = numel(R.cluster_id);
    w = collect_spike_windows(units.spike_times, units.spike_id, R.cluster_id, tr.seq_onset + win(1), tr.seq_onset + win(2));
    cnt = zeros(n, numel(w));
    for k = 1:numel(w)
        if ~isempty(w{k}); cnt(:, k) = accumarray(w{k}(:, 2), 1, [n 1]); end
    end
    ps = R.contrib_stim ./ cnt;
    ps(cnt == 0) = NaN;
    R.spikes_per_window = mean(cnt, 2);
    R.per_spike = mean(ps, 2, 'omitnan');
    R.per_spike_sem = std(ps, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(ps), 2));
    R.n_windows_active = sum(cnt > 0, 2);
    R.per_spike_p = nan(n, 1);
    for i = 1:n
        x = ps(i, ~isnan(ps(i, :)));
        if numel(x) >= 10; R.per_spike_p(i) = signrank(x); end
    end
    R.per_spike_sig = fdr_bh(R.per_spike_p, 0.05);
    R.per_spike_matrix = ps;
    pops(p) = R; %#ok<AGROW>
    fprintf('%s: per-spike contribution significant + %d, - %d of %d (FDR 0.05); per spike vs rate rho = %.2f, vs |order index - 0.5| rho = %.2f\n', ...
        R.name, sum(R.per_spike_sig & R.per_spike > 0), sum(R.per_spike_sig & R.per_spike < 0), n, ...
        corr(R.per_spike, R.rate_hz, 'type', 'Spearman', 'rows', 'complete'), ...
        corr(R.per_spike, abs(R.order_index - 0.5), 'type', 'Spearman', 'rows', 'complete'));
end

f1 = plot_rows(pops, tag);
f2 = plot_heat(pops, tag);
names = {sprintf('decoding_%s_seqlik_contrib_per_spike', stimulus_name), ...
    sprintf('decoding_%s_seqlik_contrib_per_spike_heatmap', stimulus_name)};
save_vector_figures([f1 f2], names, out_dir, cfg);
end


function sig = fdr_bh(p, q)
sig = false(size(p));
ok = find(~isnan(p));
[ps, o] = sort(p(ok));
m = numel(ps);
k = find(ps <= (1:m)' / m * q, 1, 'last');
if ~isempty(k); sig(ok(o(1:k))) = true; end
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function cmap = div_cmap()
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 0.94 0.94 0.94; 0.922 0.408 0.204], linspace(0, 1, 256));
end


function f = plot_rows(pops, tag)
np = numel(pops);
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.6 0.6 0.6];
f = figure('Color', 'w', 'Position', [30 30 1750 330 * np]);
tl = tiledlayout(f, np, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - L_{seq} contribution per spike (z units per spike, held-out sequences). ' ...
    'Total contribution = spikes per window \\times contribution per spike'], tag), 'FontSize', 11);
for p = 1:np
    R = pops(p);
    n = numel(R.cluster_id);
    pname = sprintf('%s (n = %d)', strrep(R.name, '_', ' '), n);

    % 1. per spike, ranked
    ax = nexttile(tl); hold(ax, 'on');
    [~, o] = sort(R.per_spike, 'descend', 'MissingPlacement', 'last');
    b = bar(ax, 1:n, R.per_spike(o), 1, 'FaceColor', 'flat', 'EdgeColor', 'none');
    col = repmat(grey, n, 1);
    up = R.per_spike_sig(o) & R.per_spike(o) > 0; dn = R.per_spike_sig(o) & R.per_spike(o) < 0;
    col(up, :) = repmat(orange, sum(up), 1); col(dn, :) = repmat(blue, sum(dn), 1);
    b.CData = col;
    errorbar(ax, 1:n, R.per_spike(o), R.per_spike_sem(o), 'k', 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.5);
    yline(ax, 0, 'k-'); xlim(ax, [0 n + 1]);
    xlabel(ax, 'neuron (sorted by per-spike contribution)'); ylabel(ax, 'contribution per spike (z units)');
    title(ax, sprintf('%s: per spike, %d significantly + / %d - (FDR 0.05)', pname, sum(up), sum(dn)), ...
        'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 2. decomposition: spikes x per spike = total
    ax = nexttile(tl); hold(ax, 'on');
    x = R.spikes_per_window; y = R.mean_stim;
    xl = [max(min(x(x > 0)) * 0.7, 0.05) max(x) * 1.4];
    xs = logspace(log10(xl(1)), log10(xl(2)), 50);
    yl = [min(0, min(y)) - 0.01, max(y) * 1.15];
    for e = [0.001 0.005 0.02]
        plot(ax, xs, e * xs, ':', 'Color', [0.65 0.65 0.65]);
        xe = min(xl(2), yl(2) / e) / 1.15;   % label where the line leaves the plot
        text(ax, xe, e * xe, sprintf('%.3g/spike ', e), 'FontSize', 7, 'Color', [0.45 0.45 0.45], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
    scatter(ax, x, y, 16 + 30 * R.per_spike_sig, R.per_spike, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
    ylim(ax, yl);
    lim = max(abs(R.per_spike), [], 'omitnan');
    colormap(ax, div_cmap()); clim(ax, [-lim lim]); cb = colorbar(ax); cb.Label.String = 'per spike';
    set(ax, 'XScale', 'log'); xlim(ax, xl);
    yline(ax, 0, 'k-');
    xlabel(ax, 'mean spikes per 0.8 s window'); ylabel(ax, 'mean total contribution (z units)');
    title(ax, 'decomposition: total = spikes \times per spike (dotted = equal per-spike value)', ...
        'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 3. per spike vs order index
    ax = nexttile(tl); hold(ax, 'on');
    scatter(ax, R.order_index, R.per_spike, 16 + 30 * R.per_spike_sig, log10(max(R.rate_hz, 0.1)), 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
    colormap(ax, parula); cb = colorbar(ax); cb.Label.String = 'log_{10} rate (Hz)';
    xline(ax, 0.5, 'k:'); yline(ax, 0, 'k:');
    rho = corr(R.per_spike, abs(R.order_index - 0.5), 'type', 'Spearman', 'rows', 'complete');
    xlabel(ax, 'order index (early \leftarrow 0.5 \rightarrow late)'); ylabel(ax, 'contribution per spike');
    title(ax, sprintf('per spike vs distance from the middle of the order: \\rho = %.2f', rho), ...
        'FontWeight', 'normal', 'FontSize', 8);
    style(ax);

    % 4. per spike vs order strength
    ax = nexttile(tl); hold(ax, 'on');
    scatter(ax, R.order_strength, R.per_spike, 16 + 30 * R.per_spike_sig, R.order_index, 'filled', ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
    colormap(ax, div_cmap()); cb = colorbar(ax); cb.Label.String = 'order index';
    yline(ax, 0, 'k:');
    rho = corr(R.per_spike, R.order_strength, 'type', 'Spearman', 'rows', 'complete');
    xlabel(ax, 'order strength (mean |W| with the other neurons)'); ylabel(ax, 'contribution per spike');
    title(ax, sprintf('per spike vs order strength: \\rho = %.2f', rho), 'FontWeight', 'normal', 'FontSize', 8);
    style(ax);
end
end


function f = plot_heat(pops, tag)
np = numel(pops);
f = figure('Color', 'w', 'Position', [30 30 520 * np 620]);
tl = tiledlayout(f, 1, np, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - contribution per spike in every held-out sequence (rows: neurons sorted by mean per ' ...
    'spike; columns: sequences sorted by z; white = neuron silent in that window)'], tag), 'FontSize', 11);
for p = 1:np
    R = pops(p);
    [~, o] = sort(R.per_spike, 'descend', 'MissingPlacement', 'last');
    [~, s] = sort(R.z_stim, 'descend');
    M = R.per_spike_matrix(o, s);
    lim = max(0.005, prctile(abs(M(:)), 98));
    ax = nexttile(tl);
    imagesc(ax, M, 'AlphaData', ~isnan(M), [-lim lim]); colormap(ax, div_cmap());
    cb = colorbar(ax); cb.Label.String = 'per spike (z units)';
    xlabel(ax, 'held-out sequence (sorted by z)'); ylabel(ax, 'neuron (sorted by mean per spike)');
    title(ax, sprintf('%s (n = %d)', strrep(R.name, '_', ' '), numel(R.cluster_id)), 'FontWeight', 'normal');
    set(ax, 'Color', 'w');
    style(ax);
end
end
