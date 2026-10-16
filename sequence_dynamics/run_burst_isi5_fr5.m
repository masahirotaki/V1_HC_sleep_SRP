% RUN_BURST_ISI5_FR5  Burst analysis of V1 units in the day-5 sequence tests (ABCD, A_CD, E_CD).
%
% Reproduces everything in Y:\Ellie\DATA\visual_sequence_oscillation\burst_isi5ms_fr5Hz:
%   data\                       one .mat per test day (struct r, see compute_session_burst_psth)
%   smoothed_10ms\, unsmoothed\  pooled figures per protocol (plot_burst_groups: figA-D,
%                               plot_burst_order: figE-G), each as .png, vector .pdf and .fig
%   per_animal_order\<smoothing>\  one recording-order figure per test day (plot_order_per_mouse)
%
% Parameters
%   burst          >= 2 spikes with every ISI <= 5 ms (Manuscript/bursts.py rule)
%   rate groups    narrow (cell_type 2) | broad, baseline >= 5 Hz | broad, baseline < 5 Hz,
%                  baseline = mean rate -0.6 to -0.2 s before A, averaged over the 3 recordings
%   PSTH           5 ms bins aligned to A (or E) onset; smoothing 10 ms Gaussian SD or none
%   first / last   20 trials of each recording
%
% Reads Ellie's data only; writes only to out_root. Set recompute = false to re-plot
% from the saved .mat files without touching the spike data.

recompute = true;
out_root = 'Y:\Ellie\DATA\visual_sequence_oscillation\burst_isi5ms_fr5Hz';
fr_thr = 5;                                       % Hz, broad-unit rate split
smoothing = {'smoothed_10ms', 0.010; 'unsmoothed', 0};
opt.isi_thr = 0.005;                              % s, burst ISI threshold
opt.n_jit = 1;                                    % jitter repeats (only for the unused bi_jit metric)

setup_paths();
cfg = project_config();
data_dir = fullfile(out_root, 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

%% 1. Per-session spikes -> burst labels -> PSTHs and per-trial window counts
S = day5_test_sessions();
R = cell(numel(S), 1);
for i = 1:numel(S)
    f = fullfile(data_dir, sprintf('burst_psth_%s_%s_%d.mat', S(i).subject, S(i).date, S(i).protocol));
    if ~recompute && exist(f, 'file')
        tmp = load(f);
        R{i} = tmp.r;
        continue
    end
    r = compute_session_burst_psth(S(i).subject, S(i).date, S(i).protocol, cfg, opt);
    r.week = S(i).week;
    save(f, 'r', '-v7.3');
    R{i} = r;
    fprintf('%s %s %d ms: recorded %s\n', S(i).subject, S(i).date, S(i).protocol, ...
        strjoin(r.rec_order(1, :) + " (g" + r.rec_order(3, :) + ")", ' -> '));
end
R = [R{:}];

%% 2. Figures, with and without smoothing
for k = 1:size(smoothing, 1)
    d = fullfile(out_root, smoothing{k, 1});
    for p = [150 250]
        Rp = R([R.protocol] == p);
        plot_burst_groups(Rp, d, sprintf('%dms', p), fr_thr, smoothing{k, 2});   % figA-D
        plot_burst_order(Rp, d, sprintf('%dms', p), fr_thr, smoothing{k, 2});    % figE-G
    end
    for i = 1:numel(R)
        plot_order_per_mouse(R(i), fullfile(out_root, 'per_animal_order', smoothing{k, 1}), smoothing{k, 2});
    end
end
