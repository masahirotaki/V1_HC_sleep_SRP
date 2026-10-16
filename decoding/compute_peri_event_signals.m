function S = compute_peri_event_signals(subject, date, sleep_name, cfg, opts)
% COMPUTE_PERI_EVENT_SIGNALS  Peri-event LFP / CSD (every channel of each recorded shank)
% and CA1 / V1 spiking around ripples, UP onsets and DOWN onsets of a sleep recording,
% saved once so that later analyses do not re-read the raw LFP.
%
%   S = compute_peri_event_signals('M00069', '20250929', 'Sleep_Box_2', project_config())
%
% Needs run_sleep_detection (ripples, slow oscillations, behavioural_state).
% Events (SWS only):
%   ripple_all / ripple_solo / ripple_cluster_R1   ripple onset (opts.ripple_align = 'onset', z = 1
%                                                   envelope crossing) or peak ('peak') times; solo and
%                                                   cluster after Voroslakos et al. 2026 (classify_ripple_clusters);
%                                                   S.ripple_peak_lag = peak - onset of the same events
%   up_onset / down_onset                           start / end of UP states <= opts.max_up (2 s)
% LFP: CatGT lf.bin snippets (memmapfile, no full load), uV, opts.lfp_window around each event,
%   at most opts.max_events per event type (random subset, opts.seed); each channel's mean
%   over the first 30% of the window (pre-event baseline) is subtracted. Per shank the two
%   channel columns are averaged at each depth; CSD = -d2V/dz2 on the depth-sorted mean LFP
%   (Gaussian spatial smoothing, opts.csd_sigma_um), uV/mm^2; sink < 0.
% PSD by depth: Welch PSD (Hann, opts.psd_win s) of opts.psd_n_windows random SWS windows,
%   averaged over the two columns at each depth; S.psd.shank(s).band_power for opts.psd_bands
%   (as in VR_NPX_analysis calculate_channel_PSD; plotted normalised to the shank maximum).
% Bad channels (geo(s).bad / bad_rows): zero / non-finite power, or log10 band power more than
%   opts.bad_log10 (0.3 = 2x) away from the running median of the 9 depth-nearest channels in
%   any of opts.bad_bands (<= 60 Hz, which dominate the CSD; ripple / spike bands are left out
%   because CA1 pyramidale peaks are real). They are dropped from the per-depth mean (LFP, CSD,
%   PSD); depths with no good channel are linearly interpolated across depth (geo(s).interp_depth).
% Spikes (all events): per-unit mean rate (Hz, opts.spike_bin) and per-event population
%   counts for CA1 (all units) and V1 (all units; .is_good marks good ones).
% MUA z: pooled spikes of all units per region, 1 ms bins, Gaussian sigma 5 ms, z-scored with
%   the mean / SD over all SWS (detect_mua_bursts, same as the MUA-burst detector, so z = 3 is
%   the burst threshold); peri-event mean and SEM across events.
% Output: <sleep>\detection\extracted_peri_event_signals.mat, extracted_peri_event_csd.fig/.pdf

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.max_up = 2;
def.ripple_align = 'onset';
def.lfp_window = struct('ripple', [-0.5 0.5], 'updown', [-0.5 0.5]);
def.max_events = 400;
def.csd_sigma_um = 30;
def.psd_n_windows = 150;
def.psd_win = 2;
def.psd_bands = [0.5 3; 9 17; 125 300; 300 600];
def.bad_bands = [0.5 3; 9 17; 30 60];
def.bad_log10 = 0.3;
def.spike_bin = 0.005;
def.mua_bin = 0.001;
def.mua_sigma = 0.005;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
rng(opts.seed);
paths  = get_session_paths(subject, date, sleep_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
tag = sprintf('%s %s %s', subject, date, strrep(sleep_name, '_', ' '));
B = load(fullfile(paths.output_dir, 'behavioural_state.mat')); states = B.behavioural_state;
R = load(fullfile(paths.detection_dir, 'extracted_ripple_events.mat')); rip = R.ripples;
SO = load(fullfile(paths.detection_dir, 'extracted_slow_oscillation_events.mat')); so = SO.slow_oscillations;

%% events
C = classify_ripple_clusters(rip.onset(:), rip.offset(:));
on = rip.onset(C.order); on = on(:);
pk = rip.peaktimes(C.order); pk = pk(:);
switch opts.ripple_align
    case 'onset'; rt = on;
    case 'peak';  rt = pk;
    otherwise; error('opts.ripple_align must be ''onset'' or ''peak''');
end
in_sws = @(t) any(t(:) >= states.SWS(:, 1)' & t(:) < states.SWS(:, 2)', 2);
U = so.UP(diff(so.UP, [], 2) > 0 & diff(so.UP, [], 2) <= opts.max_up, :);
sel.ripple_all = in_sws(rt);
sel.ripple_solo = in_sws(rt) & C.type == "solo";
sel.ripple_cluster_R1 = in_sws(rt) & C.position == 1;
for n = fieldnames(sel)'
    ev.(n{1}) = rt(sel.(n{1}));
    peak_lag.(n{1}) = pk(sel.(n{1})) - on(sel.(n{1}));
end
ev.up_onset = U(in_sws(U(:, 1)), 1);
ev.down_onset = U(in_sws(U(:, 2)), 2);
names = fieldnames(ev);
fprintf('%s events (SWS): %s\n', sleep_name, strjoin(cellfun(@(n) sprintf('%s %d', n, numel(ev.(n))), names, ...
    'UniformOutput', false), ', '));

%% LFP / CSD from the binary
meta = read_spikeglx_meta(paths.lf_meta_file);
n_saved = str2double(meta.nSavedChans);
fs = str2double(meta.imSampRate);
info = dir(paths.lf_bin_file);
n_samp = floor(info.bytes / (2 * n_saved));
uv = microvolts_per_bit(meta);
mm = memmapfile(paths.lf_bin_file, 'Format', {'int16', [n_saved n_samp], 'x'});
shanks = [layout.shank.id];
for s = 1:numel(shanks)
    rows = find(layout.kcoords == shanks(s));
    [dep, ~, g] = unique(layout.ycoords(rows));          % average the two columns at each depth
    geo(s).shank = shanks(s); %#ok<AGROW>
    geo(s).depth = dep;
    geo(s).rows = rows;
    geo(s).group = g;
    geo(s).bin_channels = layout.chanMap(rows);
    geo(s).CA1_depth = layout.shank(s).CA1_depth;
    geo(s).L5_depth = layout.shank(s).L5_depth;
    geo(s).L4_depth = layout.shank(s).L4_depth;
    geo(s).surface_depth = layout.shank(s).surface_depth;
end
%% PSD by depth (SWS) and bad channels
nw = round(opts.psd_win * fs);
starts = [];
for i = 1:size(states.SWS, 1)
    if diff(states.SWS(i, :)) >= opts.psd_win
        starts = [starts, states.SWS(i, 1):opts.psd_win:states.SWS(i, 2) - opts.psd_win]; %#ok<AGROW>
    end
end
if numel(starts) > opts.psd_n_windows; starts = sort(starts(randperm(numel(starts), opts.psd_n_windows))); end
a = round(starts * fs) + 1;
a = a(a >= 1 & a + nw - 1 <= n_samp);
pacc = 0;
tic;
for e = 1:numel(a)
    X = double(mm.Data.x(:, a(e) + (0:nw - 1)))';
    [pxx, fxx] = pwelch(X - mean(X, 1), hann(nw), 0, nw, fs);
    pacc = pacc + pxx;
end
pacc = pacc / max(numel(a), 1) * uv^2;                    % uV^2/Hz, freq x channel
keepf = fxx <= 700;
psd = struct('f', fxx(keepf), 'bands', opts.psd_bands, 'n_windows', numel(a));
for s = 1:numel(shanks)
    geo(s).bad = psd_bad_channels(pacc(:, geo(s).bin_channels), fxx, layout.ycoords(geo(s).rows), ...
        opts.bad_bands, opts.bad_log10);
    geo(s).bad_rows = geo(s).rows(geo(s).bad);
    geo(s).interp_depth = geo(s).depth(accumarray(geo(s).group, ~geo(s).bad) == 0);
    Pd = depth_mean(pacc(keepf, geo(s).bin_channels)', geo(s));                         % depth x freq
    psd.shank(s).power = single(Pd);
    psd.shank(s).band_power = cell2mat(arrayfun(@(b) mean(Pd(:, psd.f >= opts.psd_bands(b, 1) & ...
        psd.f <= opts.psd_bands(b, 2)), 2), 1:size(opts.psd_bands, 1), 'UniformOutput', false));
    fprintf('  shank %d: %d bad channels (rows %s), %d depths interpolated\n', geo(s).shank, numel(geo(s).bad_rows), ...
        mat2str(geo(s).bad_rows'), numel(geo(s).interp_depth));
end
fprintf('  PSD %d x %.0f s SWS windows in %.0f s\n', numel(a), opts.psd_win, toc);

%% peri-event LFP / CSD (bad channels excluded, depths without a good channel interpolated)
for k = 1:numel(names)
    nm = names{k};
    t0 = ev.(nm);
    if startsWith(nm, 'ripple'); w = opts.lfp_window.ripple; else; w = opts.lfp_window.updown; end
    if numel(t0) > opts.max_events; t0 = sort(t0(randperm(numel(t0), opts.max_events))); end
    idx = round(w(1) * fs):round(w(2) * fs);
    a = round(t0 * fs) + 1;
    a = a(a + idx(1) >= 1 & a + idx(end) <= n_samp);
    acc = zeros(n_saved, numel(idx));
    tic;
    for e = 1:numel(a)
        acc = acc + double(mm.Data.x(:, a(e) + idx));
    end
    acc = acc / max(numel(a), 1) * uv;
    % remove each channel's DC offset: subtract its mean over the pre-event baseline
    % (first 30% of the window), otherwise static channel offsets dominate the CSD
    tt = idx / fs;
    base = tt < w(1) + 0.3 * diff(w);
    acc = acc - mean(acc(:, base), 2);
    for s = 1:numel(shanks)
        V = depth_mean(acc(geo(s).bin_channels, :), geo(s));                           % depth x time
        P.(nm).shank(s).lfp = single(V);
        P.(nm).shank(s).csd = single(csd(V, geo(s).depth, opts.csd_sigma_um));
    end
    P.(nm).t = idx / fs;
    P.(nm).n_events = numel(a);
    P.(nm).event_times = t0;
    fprintf('  LFP %-18s %4d events read in %.0f s\n', nm, numel(a), toc);
end

%% spikes: per-unit mean rate and per-event population counts
for k = 1:numel(names)
    nm = names{k};
    t0 = ev.(nm);
    if startsWith(nm, 'ripple'); w = opts.lfp_window.ripple; else; w = opts.lfp_window.updown; end
    edges = w(1):opts.spike_bin:w(2) + 1e-9;
    for region = ["CA1" "V1"]
        ids = units.cluster_id(units.region == region);
        [usum, pop] = peri_counts(units.spike_times, units.spike_id, ids, t0, edges);
        Q.(nm).(region).cluster_id = ids;
        Q.(nm).(region).is_good = units.is_good(units.region == region);
        Q.(nm).(region).unit_rate = usum / max(numel(t0), 1) / opts.spike_bin;      % unit x bin (Hz)
        Q.(nm).(region).event_pop_count = single(pop);                               % event x bin
    end
    Q.(nm).t = edges(1:end - 1) + opts.spike_bin / 2;
    Q.(nm).event_times = t0;
end

%% MUA z (all units per region, SWS-referenced, same as the MUA-burst detector)
t_end = n_samp / fs;
for region = ["CA1" "V1"]
    ids = units.cluster_id(units.region == region);
    [~, mua] = detect_mua_bursts(units.spike_times(ismember(units.spike_id, ids)), t_end, ...
        struct('bin', opts.mua_bin, 'sigma', opts.mua_sigma, 'baseline', states.SWS));
    for k = 1:numel(names)
        nm = names{k};
        if startsWith(nm, 'ripple'); w = opts.lfp_window.ripple; else; w = opts.lfp_window.updown; end
        lags = round(w(1) / opts.mua_bin):round(w(2) / opts.mua_bin);
        b0 = floor(ev.(nm)(:) / opts.mua_bin) + 1;            % bin containing each event
        b0 = b0(b0 + lags(1) >= 1 & b0 + lags(end) <= numel(mua.z));
        Z = reshape(mua.z(b0 + lags), numel(b0), numel(lags));
        Q.(nm).(region).mua_z_mean = mean(Z, 1);
        Q.(nm).(region).mua_z_sem = std(Z, 0, 1) / sqrt(max(numel(b0), 1));
        Q.(nm).mua_t = lags * opts.mua_bin;
    end
end

S = struct('subject', subject, 'date', date, 'sleep', sleep_name, 'opts', opts, 'events', ev, 'geometry', geo, ...
    'lfp', P, 'spikes', Q, 'psd', psd, 'fs_lfp', fs, 'UP', U, 'clusters', C, 'ripple_peak_lag', peak_lag);
fig = plot_signals(S, tag);
if ~opts.save; return; end
safe_save(fullfile(paths.detection_dir, 'extracted_peri_event_signals.mat'), S, cfg);
save_vector_figures(fig, {'extracted_peri_event_csd'}, paths.detection_dir, cfg);
end


function [usum, pop] = peri_counts(spike_times, spike_id, ids, t0, edges)
% usum: unit x bin spike counts summed over events; pop: event x bin population counts.
% Memory-light for long recordings (per-event index ranges into the sorted spike train).
in = ismember(spike_id, ids);
[st, o] = sort(spike_times(in));
[~, col] = ismember(spike_id(in), ids);
col = col(o);
nb = numel(edges) - 1;
usum = zeros(numel(ids), nb);
pop = zeros(numel(t0), nb);
if isempty(t0) || isempty(st); return; end
[ts, ord] = sort(t0(:));
jit = (1:numel(ts))' * 1e-9;                                  % strictly increasing edges
before = @(x) cumsum(histcounts(st, [-inf; x(:) + jit]))';   % number of spikes < each x
lo = before(ts + edges(1));
hi = before(ts + edges(end));
for e = 1:numel(ts)
    r = lo(e) + 1:hi(e);
    if isempty(r); continue; end
    b = discretize(st(r) - ts(e), edges);
    ok = ~isnan(b);
    c = accumarray([col(r(ok)) b(ok)], 1, [numel(ids) nb]);
    usum = usum + c;
    pop(ord(e), :) = sum(c, 1);
end
end


function bad = psd_bad_channels(pch, f, y, bands, thr)
% pch: freq x channel PSD of one shank; y: channel depths. Bad = zero / non-finite power or a
% band power more than thr (log10) away from the running median of the 9 depth-nearest channels.
lp = cell2mat(arrayfun(@(b) log10(mean(pch(f >= bands(b, 1) & f <= bands(b, 2), :), 1))', ...
    1:size(bands, 1), 'UniformOutput', false));
dead = any(~isfinite(lp), 2);
lp(~isfinite(lp)) = NaN;
[~, o] = sort(y);
med = NaN(size(lp));
med(o, :) = movmedian(lp(o, :), 9, 1, 'omitnan');
bad = dead | any(abs(lp - med) > thr, 2);
end


function V = depth_mean(X, G)
% X: channel x time in G.bin_channels order -> depth x time mean over the good channels at each
% depth; depths with no good channel are linearly interpolated across depth (nearest at the ends)
nd = numel(G.depth);
V = NaN(nd, size(X, 2));
for d = 1:nd
    k = G.group == d & ~G.bad;
    if any(k); V(d, :) = mean(X(k, :), 1); end
end
miss = all(isnan(V), 2);
if any(miss) && sum(~miss) >= 2
    V(miss, :) = interp1(G.depth(~miss), V(~miss, :), G.depth(miss), 'linear');
    still = all(isnan(V), 2);
    if any(still); V(still, :) = interp1(G.depth(~still), V(~still, :), G.depth(still), 'nearest', 'extrap'); end
end
end


function D = csd(V, depth, sigma_um)
% -d2V/dz2 along depth (um), after Gaussian smoothing across depth; uV / mm^2
dz = median(diff(depth));
s = sigma_um / dz;
k = exp(-0.5 * ((-ceil(3 * s):ceil(3 * s)) / s).^2)'; k = k / sum(k);
Vs = conv2(V, k, 'same') ./ conv2(ones(size(V)), k, 'same');
D = zeros(size(V));
D(2:end - 1, :) = -(Vs(3:end, :) - 2 * Vs(2:end - 1, :) + Vs(1:end - 2, :)) / (dz / 1000)^2;
D([1 end], :) = NaN;
end

function c = microvolts_per_bit(meta)
range_max = str2double(meta.imAiRangeMax);
if isfield(meta, 'imMaxInt'); max_int = str2double(meta.imMaxInt); else; max_int = 512; end
if isfield(meta, 'imChan0lfGain'); gain = str2double(meta.imChan0lfGain);
elseif isfield(meta, 'imChan0apGain'); gain = str2double(meta.imChan0apGain);
else; gain = 80; end
c = range_max / max_int / gain * 1e6;
end

function f = plot_signals(S, tag)
names = fieldnames(S.lfp);
show = names(ismember(names, {'ripple_solo', 'ripple_cluster_R1', 'up_onset', 'down_onset'}));
ns = numel(S.geometry);
ncol = 5 * ns + 6;                         % per shank: PSD (2) + CSD (3); rate (3); MUA z (3)
f = figure('Color', 'w', 'Position', [30 30 140 * ncol 290 * numel(show)]);
tl = tiledlayout(f, numel(show), ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - SWS PSD by depth, peri-event CSD per shank (sink < 0, blue), CA1 / V1 rate and MUA z ' ...
    '(ripples aligned to %s)'], tag, S.opts.ripple_align), 'FontSize', 10.5);
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 0.97 0.97 0.97; 0.922 0.408 0.204], linspace(0, 1, 256));
pcol = [0.35 0.35 0.35; 0.58 0.40 0.74; 0.922 0.408 0.204; 0.10 0.62 0.55];
bands = S.psd.bands;
band_lbl = arrayfun(@(b) sprintf('\\color[rgb]{%.2f %.2f %.2f}%g-%g', pcol(b, :), bands(b, 1), bands(b, 2)), ...
    1:size(bands, 1), 'UniformOutput', false);
for k = 1:numel(show)
    nm = show{k};
    is_rip = startsWith(nm, 'ripple');
    lag = NaN;
    if is_rip && strcmp(S.opts.ripple_align, 'onset'); lag = 1000 * median(S.ripple_peak_lag.(nm)); end
    for s = 1:ns
        G = S.geometry(s);
        dz = median(diff(G.depth));
        ylims = [min(G.depth) - dz / 2, max(G.depth) + dz / 2];
        % PSD by depth (normalised band power)
        ax = nexttile(tl, [1 2]); hold(ax, 'on');
        bp = S.psd.shank(s).band_power ./ max(S.psd.shank(s).band_power, [], 1);
        for b = 1:size(bands, 1); plot(ax, bp(:, b), G.depth, '-', 'Color', pcol(b, :), 'LineWidth', 1.1); end
        if ~isempty(G.interp_depth)   % depths rebuilt by interpolation (no good channel)
            plot(ax, zeros(size(G.interp_depth)), G.interp_depth, 'k<', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
        end
        for L = {'CA1_depth', 'L5_depth', 'L4_depth'}
            if ~isnan(G.(L{1})); yline(ax, G.(L{1}), 'k:'); end
        end
        xlim(ax, [0 1.05]); ylim(ax, ylims);
        xlabel(ax, 'norm. power'); ylabel(ax, 'depth (\mum)');
        title(ax, {sprintf('PSD shank %d (Hz)', G.shank), [strjoin(band_lbl(1:2), '  ') ' ' strjoin(band_lbl(3:end), '  ')]}, ...
            'FontWeight', 'normal', 'FontSize', 6.5);
        set(ax, 'TickDir', 'out', 'FontSize', 8, 'Box', 'off');
        % CSD
        D = double(S.lfp.(nm).shank(s).csd);
        ax = nexttile(tl, [1 3]); hold(ax, 'on');
        lim = max(1, prctile(abs(D(:)), 99));
        imagesc(ax, 1000 * S.lfp.(nm).t, G.depth, D, [-lim lim]); colormap(ax, cmap); set(ax, 'YDir', 'normal');
        for L = {'CA1_depth', 'L5_depth', 'L4_depth'}
            if ~isnan(G.(L{1})); yline(ax, G.(L{1}), 'k:', strrep(erase(L{1}, '_depth'), '_', ' '), 'FontSize', 7); end
        end
        xline(ax, 0, 'k-');
        if ~isnan(lag); xline(ax, lag, 'k--'); end
        axis(ax, 'tight'); ylim(ax, ylims);
        xlabel(ax, 'ms');
        cb = colorbar(ax); cb.Label.String = 'CSD (\muV/mm^2)';
        title(ax, sprintf('%s, shank %d (n = %d)', strrep(nm, '_', ' '), G.shank, S.lfp.(nm).n_events), ...
            'FontWeight', 'normal', 'FontSize', 8);
        set(ax, 'TickDir', 'out', 'FontSize', 8);
    end
    ax = nexttile(tl, [1 3]); hold(ax, 'on');
    t = 1000 * S.spikes.(nm).t;
    rc = mean(S.spikes.(nm).CA1.unit_rate, 1); rv = mean(S.spikes.(nm).V1.unit_rate(S.spikes.(nm).V1.is_good, :), 1);
    plot(ax, t, movmean(rc, 3), '-', 'Color', [0.922 0.408 0.204], 'LineWidth', 1.6, 'DisplayName', 'CA1 (all units)');
    plot(ax, t, movmean(rv, 3), '-', 'Color', [0.165 0.471 0.839], 'LineWidth', 1.6, 'DisplayName', 'V1 (good units)');
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    if ~isnan(lag); xline(ax, lag, 'k--', 'HandleVisibility', 'off'); end
    xlabel(ax, 'ms'); ylabel(ax, 'mean rate per unit (Hz)');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s: population rate (n = %d)', strrep(nm, '_', ' '), numel(S.spikes.(nm).event_times)), ...
        'FontWeight', 'normal', 'FontSize', 8);
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
    % MUA z (SWS-referenced) with SEM band
    ax = nexttile(tl, [1 3]); hold(ax, 'on');
    t = 1000 * S.spikes.(nm).mua_t;
    regs = ["CA1" "V1"]; rcol = [0.922 0.408 0.204; 0.165 0.471 0.839];
    for r = 1:2
        m = S.spikes.(nm).(regs(r)).mua_z_mean; e = S.spikes.(nm).(regs(r)).mua_z_sem;
        fill(ax, [t fliplr(t)], [m + e fliplr(m - e)], rcol(r, :), 'EdgeColor', 'none', 'FaceAlpha', 0.25, ...
            'HandleVisibility', 'off');
        plot(ax, t, m, '-', 'Color', rcol(r, :), 'LineWidth', 1.4, 'DisplayName', regs(r) + " MUA (all units)");
    end
    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    if ~isnan(lag); xline(ax, lag, 'k--', 'HandleVisibility', 'off'); end
    xlim(ax, t([1 end]));
    xlabel(ax, 'ms'); ylabel(ax, 'MUA z (SWS mean / SD)');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    ttl = sprintf('%s: MUA z, mean \\pm SEM (1 ms, \\sigma 5 ms)', strrep(nm, '_', ' '));
    if ~isnan(lag); ttl = {ttl, sprintf('dashed: median ripple peak (+%.0f ms)', lag)}; end
    title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 8);
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
end
end
