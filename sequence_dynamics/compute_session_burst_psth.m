function R = compute_session_burst_psth(subject, date, protocol, cfg, opt)
% COMPUTE_SESSION_BURST_PSTH  Per-unit PSTHs of all / burst / single spikes for the
% day-5 test recordings (ABCD, A_CD, E_CD) of one protocol, plus burstiness metrics.
%
% Units are shared by every recording of a day (one sorting of the concatenated
% day), so a cluster_id is the same unit in all three conditions. Unit selection
% (opt.region, default V1, and Ellie's quality metrics) is taken from the ABCD recording.
%
% Spikes are labelled on each recording's full spike train with burst_labels
% (ISI <= opt.isi_thr, >= 2 spikes). Burstiness is measured in "spontaneous" grey
% periods, pooled over the 3 recordings: the ~30 s pauses between presentation
% blocks (1 s after the last sequence to 0.5 s before the next block) and the
% inter-sequence grey (opt.iti_skip s after D offset to 0.3 s before the next A).
%   fsb_spont      fraction of spikes in bursts
%   fsb_spont_jit  same after jittering every spike by +-opt.jitter
%   bi_jit         (fsb - fsb_jit) / (fsb + fsb_jit). NB: dominated by refractoriness
%                  in high-rate units (jitter creates ISIs < 2 ms), so < 0 for most cells
%   bi_spont       log2 of the spontaneous ACG at 2-6 ms over 20-40 ms (per-ms pair
%                  counts, +1 pseudo-count): ~0 Poisson, < 0 refractory/regular, > 0
%                  bursty; first-order rate invariant (Royer et al. 2012 style)
%   lv_spont       local variation of ISIs (Shinomoto 2003): 1 Poisson, >1 bursty
%   acg_spont      spontaneous autocorrelogram (pair counts), 1 ms bins up to 50 ms
%
% R.psth.<class>  [nUnit x nBin x 3] mean rate (Hz) per condition, with class
%                 all | burst (all burst spikes) | event (burst onsets) | single
% R.trial_counts.<class>{cond}  per-trial counts in R.trial_win (baseline, A, pre-B)
% R.rec_order     recording name / session_id / g number of each condition
% R.zmu / R.zsd   per-unit mean / SD of single-trial 10 ms rates over the analysis
%                 window of all conditions (z-scoring as in the manuscript figures)

if nargin < 5; opt = struct(); end
opt = set_default(opt, 'isi_thr', 0.010);
opt = set_default(opt, 'bin', 0.005);
opt = set_default(opt, 'pre', 0.6);
opt = set_default(opt, 'jitter', 0.025);
opt = set_default(opt, 'n_jit', 5);
opt = set_default(opt, 'acg_max', 0.050);
opt = set_default(opt, 'iti_skip', 0.8);   % s after D offset (post-sequence ringing)
opt = set_default(opt, 'region', "V1");    % "V1" (layers from L4 CSD sink) or "CA1" (layer = "CA1")

prefix = 'GAVNIK_';
if protocol == 250; prefix = 'GAVNIK250_'; end
conds = {'ABCD', 'A_CD', 'E_CD'};
for c = 1:numel(conds)
    paths(c) = get_session_paths(subject, date, [prefix conds{c}], cfg); %#ok<AGROW>
    trials(c) = get_abcd_trials(paths(c)); %#ok<AGROW>
end
dur = trials(1).element_duration;
post = 4 * dur + 0.8;                 % D offset + 0.8 s of the following grey
cycle = min(cellfun(@(x) median(diff(x)), {trials.seq_onset}));
assert(opt.pre + post < cycle, 'Analysis window (%.2f s) overlaps the next sequence (%.2f s)', opt.pre + post, cycle);
rel_edges = (-opt.pre:opt.bin:post)';
nB = numel(rel_edges) - 1;
% per-trial count windows (s from A onset): pre-stimulus baseline, A onset, A element,
% pre-B (ends 20 ms after B onset, before any B-evoked spikes)
trial_win = [-opt.pre, -0.2; 0.03, 0.08; 0.03, dur; dur - 0.08, dur + 0.02];
n_sub = 20;   % first / last n_sub trials of each recording get their own PSTH

layout = get_probe_layout(paths(1), cfg);
U = load_session_clusters(paths(1), layout, cfg);
opt.region = string(opt.region);
sel = find(U.region == opt.region & U.is_good);
nU = numel(sel);

%% Unit info
unit.cluster_id = U.cluster_id(sel);
unit.shank      = U.shank(sel);
unit.depth      = U.depth(sel);
unit.cell_type  = U.cell_type(sel);
unit.rel_L4     = nan(nU, 1);
unit.rel_L5     = nan(nU, 1);
unit.rel_CA1    = nan(nU, 1);
for s = layout.shank
    on = unit.shank == s.id;
    unit.rel_L4(on) = unit.depth(on) - s.L4_depth;
    unit.rel_L5(on) = unit.depth(on) - s.L5_depth;
    unit.rel_CA1(on) = unit.depth(on) - s.CA1_depth;
end
if opt.region == "V1"
    unit.layer = repmat("L4", nU, 1);           % put. L4 = within 60 um of the CSD sink
    unit.layer(unit.rel_L4 >  cfg.depth_ranges.L4(2)) = "L23";
    unit.layer(unit.rel_L4 <  cfg.depth_ranges.L4(1)) = "L56";
    unit.layer(isnan(unit.rel_L4)) = "unknown";
else
    unit.layer = repmat(opt.region, nU, 1);
end

%% Spikes per condition
classes = {'all', 'burst', 'event', 'single'};
for k = 1:numel(classes)
    psth.(classes{k}) = zeros(nU, nB, numel(conds), 'single');
    psth_sub.(classes{k}) = zeros(nU, nB, numel(conds), 2, 'single');
end
z_sum = zeros(nU, 1); z_sq = zeros(nU, 1);
n_sp = zeros(nU, 1); n_bsp = zeros(nU, 1); n_bjit = zeros(nU, 1);
lv_sum = zeros(nU, 1); lv_n = zeros(nU, 1);
n_all = zeros(nU, 1); n_ball = zeros(nU, 1);
acg_edges = 0:0.001:opt.acg_max;
acg = zeros(nU, numel(acg_edges) - 1);
fr = zeros(nU, numel(conds));
spont_dur = 0;
rng(1);
for c = 1:numel(conds)
    if c == 1
        Uc = U;
    else
        Uc = load_session_clusters(paths(c), layout, cfg);
        assert(isequal(Uc.cluster_id, U.cluster_id), 'Cluster ids differ between %s and %s', conds{1}, conds{c});
    end
    st_all = split_spikes(Uc.spike_times, Uc.spike_id, unit.cluster_id);
    T = trials(c);
    rec_dur = max(Uc.spike_times);

    % spontaneous windows: block pauses and inter-sequence grey
    k_pause = T.iti > 5;
    win = [T.seq_offset(1:end - 1) + opt.iti_skip, T.seq_onset(2:end) - 0.3];
    win(k_pause, 1) = T.seq_offset(k_pause) + 1;
    win(k_pause, 2) = T.seq_onset([false; k_pause]) - 0.5;
    win = win(diff(win, 1, 2) > 0, :);
    spont_dur = spont_dur + sum(diff(win, 1, 2));
    % map spontaneous time onto one axis with > acg_max gaps between windows
    win_shift = [0; cumsum(diff(win(1:end - 1, :), 1, 2) + 1)] - win(:, 1);

    E = T.seq_onset(:)' + rel_edges;   % [nB+1 x nTrial], non-overlapping windows
    nT = size(E, 2);
    for k = 1:numel(classes); tw.(classes{k}){c} = zeros(nU, nT, size(trial_win, 1), 'single'); end
    for u = 1:nU
        st = st_all{u};
        fr(u, c) = numel(st) / rec_dur;
        lab = burst_labels(st, opt.isi_thr);
        n_all(u) = n_all(u) + numel(st);
        n_ball(u) = n_ball(u) + nnz(lab > 0);

        % trial-aligned counts
        cls_mask = {true(size(st)), lab > 0, lab == 1, lab == 0};
        for k = 1:numel(classes)
            for w = 1:size(trial_win, 1)
                Ew = T.seq_onset(:)' + trial_win(w, :)';
                tw.(classes{k}){c}(u, 1:nT, w) = trial_counts(st(cls_mask{k}), Ew);
            end
            cnt = trial_counts(st(cls_mask{k}), E);
            psth.(classes{k})(u, :, c) = mean(cnt, 2)' / opt.bin;
            psth_sub.(classes{k})(u, :, c, 1) = mean(cnt(:, 1:n_sub), 2)' / opt.bin;
            psth_sub.(classes{k})(u, :, c, 2) = mean(cnt(:, end - n_sub + 1:end), 2)' / opt.bin;
            if k == 1
                r10 = (cnt(1:2:end - 1, :) + cnt(2:2:end, :)) / (2 * opt.bin);
                z_sum(u) = z_sum(u) + sum(r10(:));
                z_sq(u) = z_sq(u) + sum(r10(:) .^ 2);
            end
        end

        % spontaneous burstiness
        w_of = window_index(st, win);
        in_sp = w_of > 0;
        n_sp(u) = n_sp(u) + nnz(in_sp);
        n_bsp(u) = n_bsp(u) + nnz(lab(in_sp) > 0);
        for j = 1:opt.n_jit
            stj = sort(st + (2 * rand(size(st)) - 1) * opt.jitter);
            labj = burst_labels(stj, opt.isi_thr);
            n_bjit(u) = n_bjit(u) + nnz(labj(window_index(stj, win) > 0) > 0) / opt.n_jit;
        end
        sts = st(in_sp) + win_shift(w_of(in_sp));     % spontaneous spikes, windows kept apart
        isi = diff(sts);
        isi(isi > 0.9) = NaN;                          % pairs across windows
        r = (isi(1:end - 1) - isi(2:end)) ./ (isi(1:end - 1) + isi(2:end));
        lv_sum(u) = lv_sum(u) + sum(r .^ 2, 'omitnan');
        lv_n(u) = lv_n(u) + nnz(~isnan(r));
        acg(u, :) = acg(u, :) + autocorr_counts(sts, acg_edges);
    end
    n_trials(c) = nT; %#ok<AGROW>
end
n10 = sum(floor(nB / 2) * n_trials);
unit.zmu = z_sum / n10;
unit.zsd = sqrt(max(z_sq / n10 - unit.zmu .^ 2, eps));
unit.fr = fr;
unit.n_spont = n_sp;
unit.fsb_spont = n_bsp ./ max(n_sp, 1);
unit.fsb_spont_jit = n_bjit ./ max(n_sp, 1);
unit.bi_jit = (n_bsp - n_bjit) ./ max(n_bsp + n_bjit, eps);
unit.bi_jit(n_sp == 0) = NaN;
lag = acg_edges(1:end - 1) * 1000;
short = sum(acg(:, lag >= 2 & lag < 6), 2);
long = sum(acg(:, lag >= 20 & lag < 40), 2);
unit.bi_spont = log2(((short + 1) / 4) ./ ((long + 1) / 20));
unit.acg_short = short;
unit.acg_long = long;
unit.lv_spont = 3 * lv_sum ./ max(lv_n, 1);
unit.lv_spont(lv_n < 20) = NaN;
unit.fsb_all = n_ball ./ max(n_all, 1);
unit.acg_spont = acg;   % pair counts per 1 ms lag bin

R.subject = subject;
R.date = date;
R.protocol = protocol;
R.conds = conds;
R.element_duration = dur;
R.t = rel_edges(1:end - 1) + opt.bin / 2;
R.acg_t = acg_edges(1:end - 1) + 0.0005;
R.n_trials = n_trials;
R.spont_dur = spont_dur;
R.unit = unit;
R.psth = psth;
R.psth_sub = psth_sub; % .<class> [nUnit x nBin x nCond x 2]: first / last n_sub trials
R.n_sub = n_sub;
R.trial_win = trial_win;
R.trial_counts = tw;   % .<class>{cond} [nUnit x nTrial x nWin] spike counts per trial window
R.rec_order = [string({paths.stimulus}); string({paths.session_id}); string([paths.gs_number])];
R.opt = opt;
R.orientation = vertcat(trials.orientation);
layers = unique(unit.layer);
counts = arrayfun(@(l) sprintf('%s %d', l, nnz(unit.layer == l)), layers, 'UniformOutput', false);
fprintf('%s %s %d ms: %d %s units (%s), spont %.0f s\n', subject, date, protocol, nU, opt.region, ...
    strjoin(counts, ', '), spont_dur);
end


function cnt = trial_counts(st, E)
% Spike counts in the bins of every trial; E = [nEdge x nTrial] non-overlapping.
nE = size(E, 1);
c = histcounts(st, E(:))';
c(end + 1) = 0;
cnt = reshape(c, nE, []);
cnt = cnt(1:end - 1, :);
end


function w = window_index(st, win)
% Index of the [start end] row of win containing each spike (0 = none); win sorted, disjoint.
w = discretize(st, reshape(win', [], 1));
w(isnan(w) | mod(w, 2) == 0) = 0;   % even bins fall between windows
w = (w + 1) / 2;
w(w == 0.5) = 0;
end


function st = split_spikes(t, id, want)
[id, o] = sort(id(:));
t = t(o);
first = find([true; diff(id) ~= 0]);
last = [first(2:end) - 1; numel(id)];
[tf, loc] = ismember(want, id(first));
st = repmat({zeros(0, 1)}, numel(want), 1);
for k = find(tf)'
    st{k} = sort(t(first(loc(k)):last(loc(k))));
end
end


function n = autocorr_counts(st, edges)
n = zeros(1, numel(edges) - 1);
for lag = 1:numel(st)
    d = st(1 + lag:end) - st(1:end - lag);
    d = d(d < edges(end));
    if isempty(d); break; end
    n = n + histcounts(d, edges);
end
end


function opt = set_default(opt, name, value)
if ~isfield(opt, name) || isempty(opt.(name)); opt.(name) = value; end
end
