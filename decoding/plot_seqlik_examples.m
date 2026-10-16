function plot_seqlik_examples(subject, date, stimulus_name, pop_name, cfg, n_examples)
% PLOT_SEQLIK_EXAMPLES  Rasters of held-out sequences from run_abcd_sequence_likelihood,
% picked across the whole range of Lseq z (best to worst), for one population.
%
%   plot_seqlik_examples('M00069', '20250929', 'GAVNIK_ABCD_1', 'CA1_all', project_config())
%
% Each held-out sequence is drawn with neurons sorted by the order index of the matrix
% that scored it (its CV fold, trained without it), so the sort is not circular.
% Order index = mean over other neurons of P(this neuron fires after them):
% bottom = tends to fire early in the 0-0.8 s window, top = late. Dot colour = same index.
% A forward sequence looks like a lower-left to upper-right drift of spikes.
% Saved as decoding_<stimulus>_seqlik_examples_<pop> in the Decoding folder.

if nargin < 5 || isempty(cfg); cfg = project_config(); end
if nargin < 6; n_examples = 20; end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths  = get_session_paths(subject, date, stimulus_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
S = load(fullfile(out_dir, sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
R = S.results;
P = R.pop(strcmp({R.pop.name}, pop_name));
trials = R.trials;
win = R.opts.window;
ids = P.cluster_id;
n = numel(ids);
stim_w = collect_spike_windows(units.spike_times, units.spike_id, ids, trials.seq_onset + win(1), ...
    trials.seq_onset + win(2));
grey_w = collect_spike_windows(units.spike_times, units.spike_id, ids, R.grey_start, R.grey_start + diff(win));
order_idx = zeros(n, max(R.fold));
for f = 1:max(R.fold)
    M = seq_order_matrix(stim_w(R.fold ~= f), n);
    order_idx(:, f) = (sum(M.p, 2) - 0.5) / (n - 1);
end

% where the grey windows came from
gap_after = arrayfun(@(s) find(trials.seq_onset < s, 1, 'last'), R.grey_start);
gap_len = trials.seq_onset(gap_after + 1) - trials.seq_offset(gap_after);
fprintf('grey windows: %d in short inter-sequence gaps (< 5 s), %d in the %d long block pauses; %d overlap another grey window\n', ...
    sum(gap_len < 5), sum(gap_len >= 5), sum(diff(trials.seq_onset) > 5), ...
    sum(any(abs(R.grey_start - R.grey_start') < diff(win) & ~eye(numel(R.grey_start)), 2)));

z = [P.stim.z];
[~, by_z] = sort(z, 'descend');
pick = by_z(round(linspace(1, numel(z), n_examples)));
n_col = 5;
n_row = ceil(n_examples / n_col) + 1;
f = figure('Color', 'w', 'Position', [30 30 1700 290 * n_row]);
tl = tiledlayout(f, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s %s %s - %s: held-out sequences across the z range (best \\rightarrow worst); neurons sorted ' ...
    'by order index from the CV fold that did not see the sequence (bottom = early, top = late)'], subject, date, ...
    strrep(stimulus_name, '_', ' '), strrep(pop_name, '_', ' ')), 'FontSize', 10);
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 0.6 0.6 0.6; 0.922 0.408 0.204], linspace(0, 1, 256));
el = 1000 * unique([trials.element_start(2:end) trials.element_end]);
for i = 1:numel(pick)
    k = pick(i);
    raster(nexttile(tl), stim_w{k}, order_idx(:, R.fold(k)), win, el, cmap, ...
        sprintf('sequence %d: z = %.2f, pct %.3f (rank %d/%d)', k, z(k), P.stim(k).pct, find(by_z == k), numel(z)));
end
% bottom row: grey windows (best, median, worst) and the z distribution
zg = [P.grey.z];
[~, gz] = sort(zg, 'descend');
gpick = gz(round(linspace(1, numel(zg), n_col - 1)));
for i = 1:numel(gpick)
    k = gpick(i);
    raster(nexttile(tl), grey_w{k}, order_idx(:, R.fold(k)), win, [], cmap, ...
        sprintf('GREY window %d: z = %.2f', k, zg(k)));
end
ax = nexttile(tl); hold(ax, 'on');
plot(ax, sort(z, 'descend'), '-', 'Color', [0.165 0.471 0.839], 'LineWidth', 2);
plot(ax, sort(zg, 'descend'), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 2);
plot(ax, find(ismember(by_z, pick)), z(by_z(ismember(by_z, pick))), 'ko', 'MarkerSize', 4);
yline(ax, norminv(0.95), '--', 'shuffle 95%', 'FontSize', 7);
yline(ax, prctile(zg, 95), '-', 'grey 95%', 'FontSize', 7);
xlabel(ax, 'rank'); ylabel(ax, 'L_{seq} z');
legend(ax, {'held-out sequences', 'grey windows', 'shown'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, 'z of all windows, sorted', 'FontWeight', 'normal');
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);

name = sprintf('decoding_%s_seqlik_examples_%s', stimulus_name, pop_name);
save_vector_figures(f, {name}, out_dir, cfg);
end


function raster(ax, w, oi, win, el, cmap, ttl)
hold(ax, 'on');
[~, o] = sort(oi);
rank = zeros(numel(o), 1);
rank(o) = 1:numel(o);
if ~isempty(w)
    scatter(ax, 1000 * (w(:, 1) + win(1)), rank(w(:, 2)), 6, oi(w(:, 2)), 'filled');
end
colormap(ax, cmap);
clim(ax, [min(oi) max(oi)]);
if ~isempty(el); xline(ax, [0 el], 'k:'); end
xlim(ax, 1000 * win); ylim(ax, [0 numel(o) + 1]);
xlabel(ax, 'ms from window start'); ylabel(ax, 'neuron (early \rightarrow late)');
title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 7.5);
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 7);
end
