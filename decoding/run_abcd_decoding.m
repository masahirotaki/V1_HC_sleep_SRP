function results = run_abcd_decoding(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_DECODING  Decode the GAVNIK ABCD sequence from CA1 (and V1, as a positive
% control) population spiking, with 10-fold cross-validation.
%
%   results = run_abcd_decoding('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Stimulus: 4 static gratings A-B-C-D shown back to back (150 ms each for GAVNIK,
% 250 ms for GAVNIK250), always in the same order, then grey (get_abcd_trials).
% Analysis windows scale with the measured element duration (timing_defaults).
% Because the order never changes, element identity and time within the sequence
% are confounded in this recording.
%
% Analyses, per population (opts.populations):
%  1. Element identity: counts from each element's onset to its offset (opts.id_window,
%     default [0, element on-time]) -> A/B/C/D.
%     Control: four consecutive 150-ms "pseudo elements" in the grey before the sequence.
%     5-class A/B/C/D/grey: plus one grey window per sequence (A onset - 0.6 s + id_window).
%     Latency sweep: same with the 150 ms window starting -100..+450 ms after onset.
%  2. Time within the sequence: 50-ms bins from -450 to +1050 ms around A onset
%     (3 bins per element, opts.time_bin) -> decode the bin (sequence / time-cell code).
%  3. Single units: split-half sorted PSTHs (time-cell plot), time modulation tests.
%  4. Visualisation: PCA trajectory of the trial-averaged population, PLS-DA of elements.
% Decoders (fit_decoder): Poisson naive Bayes, count-conditioned multinomial Bayes
% (gain-invariant), shrinkage LDA, PCA+LDA, PLS-DA, linear SVM.
% Folds hold whole sequences; nulls shuffle labels within each sequence (decode_with_null).
% Decoder templates fit on all trials (for projecting onto sleep / ripple spiking)
% are saved in results.pop(k).templates. Cluster ids are shared across the day's
% recordings (one concatenated sort), so they index the same units in Sleep_Box_x.
%
% Outputs (decoding_* files) go to <analysis>/<subject>/<date>/<stimulus>/Decoding.

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
n_el  = numel(trials.element_names);
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));
opts = timing_defaults(opts, trials);

%% Common trial structure
% element-level samples: sequence-major (A1 B1 C1 D1 A2 ...)
el_onset = reshape(trials.element_onset', [], 1);
el_label = repmat((1:n_el)', n_seq, 1);
el_group = repelem((1:n_seq)', n_el);
pre_onset = reshape((trials.seq_onset + opts.pre_offset + trials.element_start)', [], 1);
% 5 classes: A B C D + one grey window per sequence (sequence-major, A1 B1 C1 D1 grey1 A2 ...)
cls_names = [trials.element_names, {'grey'}];
n_cls = n_el + 1;
row_el = reshape((1:n_el)' + n_cls * (0:n_seq - 1), [], 1);
row_grey = n_cls * (1:n_seq)';
cls_label = repmat((1:n_cls)', n_seq, 1);
cls_group = repelem((1:n_seq)', n_cls);
time_edges = opts.time_range(1):opts.time_bin:opts.time_range(2) + 1e-9;
time_centres = time_edges(1:end - 1) + opts.time_bin / 2;
n_tb = numel(time_centres);
psth_edges = opts.psth_range(1):opts.psth_bin:opts.psth_range(2) + 1e-9;
psth_t = psth_edges(1:end - 1) + opts.psth_bin / 2;

figs = gobjects(0);
names = {};
[figs(end + 1), names{end + 1}] = plot_stimulus_structure(trials, tag);

for p = 1:numel(opts.populations)
    P = opts.populations(p);
    sel = units.region == P.region & (units.is_good | ~P.good_only);
    ids = units.cluster_id(sel);
    fprintf('\n=== %s: %d units\n', P.name, numel(ids));

    % sequence-aligned counts: trial x neuron x time bin
    counts_time = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset, time_edges);
    active = squeeze(sum(counts_time, [1 3])) >= opts.min_spikes;
    ids = ids(active);
    counts_time = counts_time(:, active, :);
    fprintf('%d units with >= %d spikes in the sequence windows\n', numel(ids), opts.min_spikes);
    counts_el  = bin_trial_spikes(units.spike_times, units.spike_id, ids, el_onset, opts.id_window);
    counts_pre = bin_trial_spikes(units.spike_times, units.spike_id, ids, pre_onset, [0 opts.id_window(2) - opts.id_window(1)]);
    counts_psth = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset, psth_edges);
    counts_grey = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset + opts.grey_offset, opts.id_window);
    counts_cls = zeros(n_cls * n_seq, numel(ids));
    counts_cls(row_el, :) = counts_el;
    counts_cls(row_grey, :) = counts_grey;

    pop = struct();
    pop.name = P.name;
    pop.region = P.region;
    pop.good_only = P.good_only;
    pop.cluster_id = ids;
    [~, row] = ismember(ids, units.cluster_id);
    pop.shank = units.shank(row);
    pop.depth = units.depth(row);
    pop.cell_type = units.cell_type(row);
    pop.counts_time = counts_time;   % sequence x neuron x time bin (opts.time_bin, time_centres)
    pop.counts_element = reshape(counts_el, n_el, n_seq, []); % element x sequence x neuron
    pop.counts_element = permute(pop.counts_element, [2 3 1]); % sequence x neuron x element
    pop.counts_pre = permute(reshape(counts_pre, n_el, n_seq, []), [2 3 1]);
    pop.counts_grey = counts_grey;   % sequence x neuron, grey window (grey_offset + id_window)
    pop.class_names = cls_names;

    %% 1. Element identity
    fprintf('-- element identity, window %s s\n', mat2str(opts.id_window));
    d = struct('n_folds', opts.n_folds, 'n_perm', opts.n_perm, 'seed', opts.seed);
    pop.identity = decode_with_null(counts_el, el_label, el_group, opts.methods, d);
    fprintf('-- control: pseudo elements in the grey (%.2f s before A)\n', -opts.pre_offset);
    pop.identity_pre = decode_with_null(counts_pre, el_label, el_group, opts.methods, d);
    fprintf('-- A/B/C/D/grey (grey window %.2f s before A)\n', -opts.grey_offset);
    pop.identity_grey = decode_with_null(counts_cls, cls_label, cls_group, opts.methods, d);

    sweep_methods = {'bayes', 'multinomial', 'lda'};
    fprintf('-- latency sweep (%s)\n', strjoin(sweep_methods, ', '));
    nm = numel(sweep_methods);
    sweep = struct('latency', opts.latencies, 'acc', nan(nm, numel(opts.latencies)), ...
        'null95', nan(nm, numel(opts.latencies)), 'p', nan(nm, numel(opts.latencies)), 'methods', {sweep_methods});
    d.n_perm = opts.n_perm_latency;
    for i = 1:numel(opts.latencies)
        c = bin_trial_spikes(units.spike_times, units.spike_id, ids, el_onset, opts.latencies(i) + [0 trials.element_duration]);
        r = evalc_quiet(@() decode_with_null(c, el_label, el_group, sweep.methods, d));
        sweep.acc(:, i) = [r.acc]';
        sweep.null95(:, i) = arrayfun(@(x) prctile(x.null_acc, 95), r)';
        sweep.p(:, i) = [r.p]';
    end
    pop.latency_sweep = sweep;

    %% 2. Time within the sequence
    fprintf('-- time bin within sequence (%d x %.0f ms bins)\n', n_tb, 1000 * opts.time_bin);
    Xt = reshape(permute(counts_time, [3 1 2]), n_tb * n_seq, []);
    yt = repmat((1:n_tb)', n_seq, 1);
    gt = repelem((1:n_seq)', n_tb);
    d = struct('n_folds', opts.n_folds, 'n_perm', opts.n_perm_time, 'seed', opts.seed, 'bin_width', opts.time_bin);
    pop.time = decode_with_null(Xt, yt, gt, opts.methods, d);
    in_stim = time_centres > 0 & time_centres < trials.seq_duration;
    for m = 1:numel(pop.time)
        ok = in_stim(yt)';
        pop.time(m).acc_stim = mean(pop.time(m).pred(ok) == yt(ok));
        pop.time(m).acc_within1_stim = mean(abs(pop.time(m).pred(ok) - yt(ok)) <= 1);
        pop.time(m).err_per_bin = accumarray(yt, abs(pop.time(m).pred - yt), [n_tb 1], @mean) * opts.time_bin;
    end

    %% 3. Single units
    pop.units = unit_time_stats(counts_time, counts_psth, time_centres, psth_t, in_stim, opts);

    %% 4. PCA trajectory and PLS-DA
    pop.pca = population_pca(counts_psth, psth_t, opts);
    pop.pls = element_pls(counts_el, el_label);

    %% Templates on all trials (for sleep projection)
    tmpl.element_window = opts.id_window;
    tmpl.element_rate_hz = mean(pop.counts_element, 1) / diff(opts.id_window); % 1 x neuron x element
    tmpl.element_rate_hz = squeeze(tmpl.element_rate_hz)';                     % element x neuron
    tmpl.time_centres = time_centres;
    tmpl.time_rate_hz = squeeze(mean(counts_time, 1))' / opts.time_bin;        % time bin x neuron
    tmpl.element_multinomial = fit_decoder(counts_el, el_label, n_el, 'multinomial'); % gain-invariant
    tmpl.time_multinomial = fit_decoder(Xt, yt, n_tb, 'multinomial');
    tmpl.element_lda = fit_decoder(counts_el, el_label, n_el, 'lda');
    tmpl.class_names = cls_names;                 % A B C D grey ("no stimulus" class)
    tmpl.grey_multinomial = fit_decoder(counts_cls, cls_label, n_cls, 'multinomial');
    tmpl.grey_lda = fit_decoder(counts_cls, cls_label, n_cls, 'lda');
    tmpl.time_lda = fit_decoder(Xt, yt, n_tb, 'lda');
    tmpl.note = ['Bayes: lambda = rate_hz * bin width of the test data; ' ...
        'LDA models expect counts in the training window width (apply_decoder).'];
    pop.templates = tmpl;

    pops(p) = pop; %#ok<AGROW>
end

%% Figures
[figs(end + 1), names{end + 1}] = plot_identity(pops, trials, opts, tag);
[figs(end + 1), names{end + 1}] = plot_time_decoding(pops, time_centres, trials, opts, tag);
[figs(end + 1), names{end + 1}] = plot_time_cells(pops, psth_t, trials, tag);
[figs(end + 1), names{end + 1}] = plot_pca_pls(pops, trials, tag);

%% Save
results.subject = subject;
results.date = date;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.trials = trials;
results.time_centres = time_centres;
results.psth_t = psth_t;
results.pop = pops;
results.summary = summary_table(pops);
disp(results.summary);
if ~opts.save; return; end

base =sprintf('decoding_%s', stimulus_name);
safe_save(fullfile(out_dir, [base '_results.mat']), struct('results', results), cfg);
csv_file = fullfile(out_dir, [base '_summary.csv']);
check_output_path(csv_file, cfg);
check_overwrite(csv_file, cfg);
writetable(results.summary, csv_file);
for i = 1:numel(figs)
    names{i} = sprintf('%s_%s', base, names{i});
end
save_vector_figures(figs, names, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function opts = fill_defaults(opts)
def.populations = struct('name', {'CA1_good', 'CA1_all', 'V1_good'}, 'region', {"CA1", "CA1", "V1"}, ...
    'good_only', {true, false, true});
def.methods = {'bayes', 'multinomial', 'lda', 'pca_lda', 'pls', 'svm'};
def.n_folds = 10;
def.n_perm = 200;                        % element identity nulls
def.n_perm_time = [200 200 50 50 50 20]; % time-bin nulls, per method (svm is slow with 20 classes)
def.n_perm_latency = 100;
% Windows left empty are set from the measured element duration (timing_defaults), so
% GAVNIK (150 ms) and GAVNIK250 (250 ms) sessions both work:
def.id_window = [];                  % [0, element on-time]: the element's actual onset -> offset
                                     % (until 2026-10-15: [0.03, 0.18] for V1 latency)
def.pre_offset = [];                 % -(sequence + 0.3): pseudo elements end 0.3 s before A (in the grey)
def.grey_offset = -0.6;              % grey class: id_window after (A onset - 0.6 s)
def.latencies = -0.1:0.025:0.45;
def.time_bin = 0.05;                 % s; 3 bins per 150-ms element (was 0.075 until 2026-10-13)
def.time_range = [];                 % [-0.45, sequence + 0.45], whole bins
def.psth_bin = 0.01;
def.psth_range = [];                 % [-0.5, sequence + 0.5]
def.psth_sigma = 0.025;              % s, Gaussian smoothing of PSTHs
def.min_spikes = 20;                 % per unit, summed over all sequence windows
def.seed = 1;
def.save = true;                     % false = compute and return without writing
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
end


function opts = timing_defaults(opts, trials)
% fill windows that depend on the stimulus timing (GAVNIK 150 ms, GAVNIK250 250 ms back to
% back; TRAIN 150 ms elements every 300 ms)
d = round(trials.element_duration, 2);      % nominal element on-time (0.15 / 0.25 s)
seq = round(trials.seq_duration, 2);        % nominal A onset -> D offset (0.6 / 1.0 / 1.05 s)
if isempty(opts.id_window); opts.id_window = [0 d]; end
if isempty(opts.pre_offset); opts.pre_offset = -(seq + 0.3); end
if isempty(opts.time_range)
    opts.time_range = [-0.45, -0.45 + opts.time_bin * ceil((seq + 0.9) / opts.time_bin)];
end
if isempty(opts.psth_range); opts.psth_range = [-0.5, seq + 0.5]; end
end


function r = evalc_quiet(fun)
% run FUN without its console output
[~, r] = evalc('fun()');
end


function U = unit_time_stats(counts_time, counts_psth, time_centres, psth_t, in_stim, opts)
% Per unit: stimulus modulation (0-600 ms vs 450 ms before, sign-rank on trial rates),
% time modulation within the sequence (Kruskal-Wallis across the 8 stimulus bins),
% split-half peak times of the smoothed PSTH (odd vs even sequences).
n = size(counts_time, 2);
pre = time_centres < 0;
rate_pre  = sum(counts_time(:, :, pre), 3) / (sum(pre) * opts.time_bin);
rate_stim = sum(counts_time(:, :, in_stim), 3) / (sum(in_stim) * opts.time_bin);
U.rate_pre = mean(rate_pre, 1)';
U.rate_stim = mean(rate_stim, 1)';
U.p_stim = nan(n, 1);
U.p_time = nan(n, 1);
stim_bins = find(in_stim);
for k = 1:n
    U.p_stim(k) = signrank(rate_stim(:, k), rate_pre(:, k));
    c = squeeze(counts_time(:, k, stim_bins));
    U.p_time(k) = kruskalwallis(c(:), repelem((1:numel(stim_bins))', size(c, 1)), 'off');
end
sm = smooth_psth(counts_psth, opts);              % trial x neuron x time (Hz)
odd = 1:2:size(sm, 1);
even = 2:2:size(sm, 1);
U.psth_mean = squeeze(mean(sm, 1));                % neuron x time
U.psth_odd  = squeeze(mean(sm(odd, :, :), 1));
U.psth_even = squeeze(mean(sm(even, :, :), 1));
[~, i_odd]  = max(U.psth_odd, [], 2);
[~, i_even] = max(U.psth_even, [], 2);
U.peak_odd  = psth_t(i_odd)';
U.peak_even = psth_t(i_even)';
U.split_half_r = arrayfun(@(k) corr(U.psth_odd(k, :)', U.psth_even(k, :)'), (1:n)');
U.alpha = 0.01;
fprintf('   units: stimulus-modulated %d/%d, time-modulated within the sequence %d/%d (p < %.2f)\n', ...
    sum(U.p_stim < U.alpha), n, sum(U.p_time < U.alpha), n, U.alpha);
end


function sm = smooth_psth(counts, opts)
sig = opts.psth_sigma / opts.psth_bin;
k = exp(-0.5 * ((-ceil(3 * sig):ceil(3 * sig)) / sig).^2);
k = reshape(k / sum(k), 1, 1, []);
% normalised convolution: no drop-off at the window edges
w = convn(ones(1, 1, size(counts, 3)), k, 'same');
sm = convn(counts, k, 'same') ./ w / opts.psth_bin;
end


function S = population_pca(counts_psth, psth_t, opts)
% PCA of the trial-averaged, soft-normalised population PSTH (rate / (range + 5 Hz)).
R = squeeze(mean(smooth_psth(counts_psth, opts), 1));     % neuron x time
R = R ./ (max(R, [], 2) - min(R, [], 2) + 5);
R = R - mean(R, 2);
[coeff, score, latent] = pca(R');
S.t = psth_t;
S.score = score(:, 1:3);
S.coeff = coeff(:, 1:3);
S.explained = 100 * latent(1:3) / sum(latent);
end


function S = element_pls(counts_el, el_label)
% PLS-DA on z-scored element-window counts (all trials; for display only - the
% cross-validated accuracy is in .identity).
Z = zscore(counts_el);
Z(:, any(isnan(Z), 1)) = 0;
Y = double(el_label == 1:max(el_label));
[~, ~, xs, ~, ~, pctvar] = plsregress(Z, Y, 3);
S.score = xs;
S.label = el_label;
S.pct_var_x = 100 * pctvar(1, :);
S.pct_var_y = 100 * pctvar(2, :);
end


function T = summary_table(pops)
rows = {};
for p = 1:numel(pops)
    sets = {'identity', 'identity_pre', 'identity_grey', 'time'};
    for s = 1:numel(sets)
        R = pops(p).(sets{s});
        for m = 1:numel(R)
            mae = NaN; acc_stim = NaN;
            if strcmp(sets{s}, 'time'); mae = R(m).mae; acc_stim = R(m).acc_stim; end
            rows(end + 1, :) = {pops(p).name, numel(pops(p).cluster_id), sets{s}, R(m).method, R(m).acc, ...
                R(m).chance, prctile(R(m).null_acc, 95), R(m).p, numel(R(m).null_acc), acc_stim, mae}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'population', 'n_units', 'analysis', 'decoder', 'accuracy', ...
    'chance', 'null_95pct', 'p_perm', 'n_perm', 'accuracy_stim_bins', 'mean_abs_error_s'});
end


%% ------------------------------------------------------------------------
%% Plotting
function c = el_colors()
% categorical slots blue, orange, aqua, violet (+ marker shape as second cue)
c = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655];
end

function cmap = seq_cmap()
% single-hue sequential ramp, white -> blue -> dark blue
cmap = interp1([0 0.6 1], [1 1 1; 0.165 0.471 0.839; 0.04 0.15 0.35], linspace(0, 1, 256));
end

function shade_elements(ax, trials, yl)
c = el_colors();
for e = 1:numel(trials.element_names)
    x0 = trials.element_start(e);
    patch(ax, [x0 x0 + trials.element_duration x0 + trials.element_duration x0], yl([1 1 2 2]), c(e, :), ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax, x0 + trials.element_duration / 2, yl(2), trials.element_names{e}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end


function [f, name] = plot_stimulus_structure(trials, tag)
name = 'stimulus_structure';
f = figure('Color', 'w', 'Position', [50 50 1100 330]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
title(tl, sprintf('%s - stimulus timing (Task\\_info.stim\\_onset = element onsets)', tag), 'FontSize', 10);
c = el_colors();
ax = nexttile(tl); hold(ax, 'on');
t0 = trials.seq_onset(1);
for e = 1:numel(trials.element_names)
    x = 1000 * ([trials.element_onset(1, e) trials.element_offset(1, e)] - t0);
    patch(ax, x([1 2 2 1]), [0 0 1 1], c(e, :), 'EdgeColor', 'w', 'LineWidth', 1);
    text(ax, mean(x), 0.5, sprintf('%s\n%d%s', trials.element_names{e}, trials.orientation(e), char(176)), ...
        'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
end
text(ax, 1000 * (trials.seq_offset(1) - t0) + 20, 0.5, sprintf('grey\n%.2f s', median(trials.iti)), 'FontSize', 8, ...
    'Color', [0.3 0.3 0.3]);
xlabel(ax, 'time from A onset (ms)'); yticks(ax, []); xlim(ax, [-100 900]);
title(ax, 'one sequence (orientation of each element)', 'FontWeight', 'normal');
style(ax);
ax = nexttile(tl);
histogram(ax, 1000 * trials.soa(:), 'BinWidth', 1, 'FaceColor', c(1, :), 'EdgeColor', 'w');
xlabel(ax, 'element onset-to-onset (ms)'); ylabel(ax, 'count');
title(ax, sprintf('element duration, median %.1f ms', 1000 * trials.element_duration), 'FontWeight', 'normal');
style(ax);
ax = nexttile(tl);
iti = trials.iti(trials.iti < 5);
histogram(ax, iti, 'BinWidth', 0.005, 'FaceColor', c(1, :), 'EdgeColor', 'w');
xlabel(ax, 'grey between sequences, D offset to next A (s)'); ylabel(ax, 'count');
title(ax, sprintf('%d sequences in %d blocks; grey %.2f s', numel(trials.seq_onset), max(trials.block), median(iti)), ...
    'FontWeight', 'normal');
style(ax);
end


function [f, name] = plot_identity(pops, trials, opts, tag)
name = 'element_identity';
np = numel(pops);
f = figure('Color', 'w', 'Position', [50 50 1250 300 * np]);
tl = tiledlayout(f, np, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - element identity (A/B/C/D), %d-fold CV, window %d-%d ms after element onset', tag, ...
    opts.n_folds, round(1000 * opts.id_window)), 'FontSize', 10);
grey = [0.6 0.6 0.6];
for p = 1:np
    P = pops(p);
    ax = nexttile(tl); hold(ax, 'on');
    acc = [[P.identity.acc]' [P.identity_pre.acc]'];
    b = bar(ax, 100 * acc, 'BarWidth', 0.8, 'EdgeColor', 'none');
    b(1).FaceColor = [0.165 0.471 0.839]; b(2).FaceColor = grey;
    n95 = [arrayfun(@(r) prctile(r.null_acc, 95), P.identity)' arrayfun(@(r) prctile(r.null_acc, 95), P.identity_pre)'];
    for k = 1:2
        plot(ax, b(k).XEndPoints, 100 * n95(:, k), '_', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
    for m = 1:numel(P.identity)
        if P.identity(m).p < 0.05
            text(ax, b(1).XEndPoints(m), 100 * acc(m, 1) + 3, '*', 'HorizontalAlignment', 'center', 'FontSize', 11);
        end
    end
    yline(ax, 25, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    set(ax, 'XTick', 1:numel(P.identity), 'XTickLabel', strrep({P.identity.method}, '_', '+'));
    ylim(ax, [0 118]); yticks(ax, 0:25:100); ylabel(ax, 'accuracy (%), dotted = chance');
    legend(ax, {'stimulus', 'grey control'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s (n = %d units): - = null 95%%, * p < 0.05', strrep(P.name, '_', ' '), numel(P.cluster_id)), ...
        'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl);
    G = P.identity_grey;
    [~, best] = max([G.acc]);
    nc = numel(P.class_names);
    imagesc(ax, 100 * G(best).confusion, [0 100]); colormap(ax, seq_cmap());
    axis(ax, 'square');
    set(ax, 'XTick', 1:nc, 'XTickLabel', P.class_names, 'YTick', 1:nc, 'YTickLabel', P.class_names);
    for i = 1:nc
        for j = 1:nc
            v = 100 * G(best).confusion(i, j);
            text(ax, j, i, sprintf('%.0f', v), 'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', [1 1 1] * (v > 55));
        end
    end
    xlabel(ax, 'decoded'); ylabel(ax, 'true');
    title(ax, sprintf('A/B/C/D/grey confusion (%%), %s: %.0f%% (chance 20%%)', ...
        strrep(G(best).method, '_', '+'), 100 * G(best).acc), 'FontWeight', 'normal');
    cb = colorbar(ax); cb.Label.String = '% of true'; style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    S = P.latency_sweep;
    lc = [0.165 0.471 0.839; 0.106 0.686 0.478; 0.922 0.408 0.204];
    lat_ms = 1000 * S.latency;
    for m = 1:numel(S.methods)
        plot(ax, lat_ms, 100 * S.null95(m, :), '--', 'Color', lc(m, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, lat_ms, 100 * S.acc(m, :), '-o', 'Color', lc(m, :), 'MarkerFaceColor', lc(m, :), ...
            'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', S.methods{m});
    end
    yline(ax, 25, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, 0, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    xlabel(ax, 'window start re element onset (ms), 150 ms window'); ylabel(ax, 'accuracy (%)');
    ylim(ax, [0 100]);
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, 'latency sweep (dashed = null 95%)', 'FontWeight', 'normal');
    style(ax);
end
end


function [f, name] = plot_time_decoding(pops, tc, trials, opts, tag)
name = 'time_in_sequence';
np = numel(pops);
f = figure('Color', 'w', 'Position', [50 50 1250 300 * np]);
tl = tiledlayout(f, np, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - decoding time within the sequence (%d ms bins, %d-fold CV); stimulus 0-%d ms', tag, ...
    round(1000 * opts.time_bin), opts.n_folds, round(1000 * trials.seq_duration)), ...
    'FontSize', 10);
n_tb = numel(tc);
% expected |error| of uniformly random guesses for each true bin
chance_err = arrayfun(@(b) mean(abs((1:n_tb) - b)), 1:n_tb) * opts.time_bin;
lc = [0.165 0.471 0.839; 0.922 0.408 0.204];
for p = 1:np
    P = pops(p);
    ax = nexttile(tl);
    imagesc(ax, 1000 * tc, 1000 * tc, 100 * P.time(1).confusion, [0 max(30, 100 * max(P.time(1).confusion(:)))]);
    colormap(ax, seq_cmap()); axis(ax, 'square'); set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    plot(ax, 1000 * tc([1 end]), 1000 * tc([1 end]), ':', 'Color', [0.5 0.5 0.5]);
    for x = [0 1000 * trials.seq_duration]
        xline(ax, x, '-', 'Color', [0.5 0.5 0.5]); yline(ax, x, '-', 'Color', [0.5 0.5 0.5]);
    end
    xlabel(ax, 'decoded time (ms)'); ylabel(ax, 'true time (ms)');
    cb = colorbar(ax); cb.Label.String = '% of true';
    title(ax, sprintf('%s (n = %d), Bayes', strrep(P.name, '_', ' '), numel(P.cluster_id)), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    yl = [0 0.65];
    shade_elements(ax, trials, yl);
    plot(ax, 1000 * tc, chance_err, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'chance');
    for m = 1:2
        plot(ax, 1000 * tc, P.time(m).err_per_bin, '-o', 'Color', lc(m, :), 'MarkerFaceColor', lc(m, :), ...
            'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', P.time(m).method);
    end
    % make the shading span ms
    for h = findobj(ax, 'Type', 'patch')'
        h.XData = 1000 * h.XData;
    end
    for h = findobj(ax, 'Type', 'text')'
        h.Position(1) = 1000 * h.Position(1);
    end
    ylim(ax, yl); xlim(ax, 1000 * tc([1 end]) + [-40 40]);
    xlabel(ax, 'true time from A onset (ms)'); ylabel(ax, 'mean |decoded - true| (s)');
    legend(ax, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 7);
    title(ax, 'decoding error per time bin', 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    acc = 100 * [P.time.acc];
    b = bar(ax, acc, 'BarWidth', 0.6, 'FaceColor', [0.165 0.471 0.839], 'EdgeColor', 'none');
    n95 = arrayfun(@(r) prctile(r.null_acc, 95), P.time);
    plot(ax, b.XEndPoints, 100 * n95, '_', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 1.5);
    plot(ax, b.XEndPoints, 100 * [P.time.acc_stim], 'd', 'Color', [0.922 0.408 0.204], ...
        'MarkerFaceColor', [0.922 0.408 0.204], 'MarkerSize', 6);
    for m = 1:numel(P.time)
        text(ax, b.XEndPoints(m), acc(m) + 2, sprintf('%.0f', acc(m)), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    yline(ax, 100 / n_tb, ':', 'chance', 'Color', [0.3 0.3 0.3], 'FontSize', 7);
    set(ax, 'XTick', 1:numel(P.time), 'XTickLabel', strrep({P.time.method}, '_', '+'));
    ylabel(ax, 'exact-bin accuracy (%)'); ylim(ax, [0 max(40, 1.2 * max(acc))]);
    legend(ax, {'all bins', 'null 95%', 'stimulus bins only'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, 'accuracy per decoder', 'FontWeight', 'normal');
    style(ax);
end
end


function [f, name] = plot_time_cells(pops, t, trials, tag)
name = 'time_cells';
np = numel(pops);
f = figure('Color', 'w', 'Position', [50 50 420 * np 900]);
tl = tiledlayout(f, 3, np, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - single units aligned to A onset (sorted by odd-sequence peak, shown: even sequences)', tag), ...
    'FontSize', 10);
for p = 1:np
    U = pops(p).units;
    [~, order] = sort(U.peak_odd);
    E = U.psth_even(order, :);
    E = (E - min(E, [], 2)) ./ max(max(E, [], 2) - min(E, [], 2), eps);
    ax = nexttile(tl, p);
    imagesc(ax, 1000 * t, 1:size(E, 1), E, [0 1]); colormap(ax, seq_cmap());
    xline(ax, 0, 'k-'); xline(ax, 1000 * unique([trials.element_start(2:end) trials.element_end]), 'k:');
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'unit (sorted)');
    title(ax, sprintf('%s: time-modulated %d/%d', strrep(pops(p).name, '_', ' '), sum(U.p_time < U.alpha), ...
        numel(U.p_time)), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl, np + p); hold(ax, 'on');
    pop_rate = mean(U.psth_mean, 1);
    yl = [0 1.15 * max(pop_rate)];
    shade_elements(ax, trials, yl);
    for h = findobj(ax, 'Type', 'patch')'; h.XData = 1000 * h.XData; end
    for h = findobj(ax, 'Type', 'text')'; h.Position(1) = 1000 * h.Position(1); end
    plot(ax, 1000 * t, pop_rate, 'k', 'LineWidth', 2);
    ylim(ax, yl); xlim(ax, 1000 * t([1 end]));
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'mean rate per unit (Hz)');
    title(ax, sprintf('stimulus-modulated %d/%d (sequence vs before A)', sum(U.p_stim < U.alpha), numel(U.p_stim)), ...
        'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl, 2 * np + p); hold(ax, 'on');
    sig = U.p_time < U.alpha;
    scatter(ax, 1000 * U.peak_odd(~sig), 1000 * U.peak_even(~sig), 14, [0.7 0.7 0.7], 'filled');
    scatter(ax, 1000 * U.peak_odd(sig), 1000 * U.peak_even(sig), 18, [0.165 0.471 0.839], 'filled');
    plot(ax, 1000 * t([1 end]), 1000 * t([1 end]), 'k:');
    r = NaN;
    if sum(sig) > 3; r = corr(U.peak_odd(sig), U.peak_even(sig), 'Type', 'Spearman'); end
    axis(ax, 'square'); xlim(ax, 1000 * t([1 end])); ylim(ax, 1000 * t([1 end]));
    xlabel(ax, 'peak time, odd sequences (ms)'); ylabel(ax, 'peak time, even sequences (ms)');
    legend(ax, {'not time-modulated', 'time-modulated'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('split-half peak consistency, \\rho = %.2f (modulated)', r), 'FontWeight', 'normal');
    style(ax);
end
end


function [f, name] = plot_pca_pls(pops, trials, tag)
name = 'pca_pls';
np = numel(pops);
f = figure('Color', 'w', 'Position', [50 50 1250 300 * np]);
tl = tiledlayout(f, np, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - population structure (PCA of trial-averaged PSTH; PLS-DA of element windows, not cross-validated)', ...
    tag), 'FontSize', 10);
c = el_colors();
mk = {'o', 's', '^', 'd'};
pcc = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478];
for p = 1:np
    S = pops(p).pca;
    ax = nexttile(tl); hold(ax, 'on');
    yl = [min(S.score(:)) max(S.score(:))] * 1.15;
    shade_elements(ax, trials, yl);
    for h = findobj(ax, 'Type', 'patch')'; h.XData = 1000 * h.XData; end
    for h = findobj(ax, 'Type', 'text')'; h.Position(1) = 1000 * h.Position(1); end
    for k = 1:3
        plot(ax, 1000 * S.t, S.score(:, k), 'Color', pcc(k, :), 'LineWidth', 2, ...
            'DisplayName', sprintf('PC%d (%.0f%%)', k, S.explained(k)));
    end
    ylim(ax, yl); xlim(ax, 1000 * S.t([1 end]));
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'PC score');
    legend(ax, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s (n = %d): PCs over time', strrep(pops(p).name, '_', ' '), numel(pops(p).cluster_id)), ...
        'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, S.score(:, 1), S.score(:, 2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1, 'DisplayName', 'grey');
    for e = 1:numel(trials.element_names)
        in_e = S.t >= trials.element_start(e) & S.t < trials.element_end(e);
        plot(ax, S.score(in_e, 1), S.score(in_e, 2), '-', 'Color', c(e, :), 'LineWidth', 2.5, ...
            'DisplayName', trials.element_names{e});
    end
    [~, i0] = min(abs(S.t));
    plot(ax, S.score(i0, 1), S.score(i0, 2), 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'A onset');
    plot(ax, S.score(1, 1), S.score(1, 2), 'ko', 'MarkerSize', 6, 'DisplayName', sprintf('%d ms', round(1000 * S.t(1))));
    xlabel(ax, 'PC1'); ylabel(ax, 'PC2'); axis(ax, 'square');
    legend(ax, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 7);
    title(ax, 'trajectory (grey = outside the sequence)', 'FontWeight', 'normal');
    style(ax);

    L = pops(p).pls;
    ax = nexttile(tl); hold(ax, 'on');
    for e = 1:numel(trials.element_names)
        in_e = L.label == e;
        scatter(ax, L.score(in_e, 1), L.score(in_e, 2), 12, c(e, :), mk{e}, 'filled', 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', trials.element_names{e});
    end
    xlabel(ax, sprintf('PLS LV1 (%.0f%% of Y)', L.pct_var_y(1))); ylabel(ax, sprintf('PLS LV2 (%.0f%% of Y)', L.pct_var_y(2)));
    axis(ax, 'square');
    legend(ax, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 7);
    title(ax, 'PLS-DA scores, one dot per element presentation', 'FontWeight', 'normal');
    style(ax);
end
end
