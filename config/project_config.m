function cfg = project_config()
% PROJECT_CONFIG  Paths and analysis parameters for the V1-HC sleep SRP project.
%
% Outputs are written into each recording's existing analysis folder on the server,
% SUBJECTS/<subj>/analysis/<date>/<stimulus>/ (+ detection/). Everything else on
% Y:\ is read-only. check_output_path / check_overwrite enforce this: only those
% folders are writable, and an existing file is only overwritten if its name is
% one of this pipeline's outputs (cfg.output.own_files).

%% Paths
cfg.data_root     = 'Y:\Ellie\DATA';
cfg.subjects_root = fullfile(cfg.data_root, 'SUBJECTS');
cfg.stim_table    = fullfile(cfg.data_root, 'Ellie_session_stimuli_table.xlsx');
cfg.output_root   = cfg.subjects_root;
% Existing files with these names may be overwritten. In the stimulus folder this
% replaces Ellie's earlier behavioural_state.mat (user decision 2026-10-13).
cfg.output.own_files.stimulus  = {'behavioural_state.*', 'sleep_stages.*', 'extracted_LFP.mat'};
cfg.output.own_files.detection = {'extracted_*', 'slow_oscillation_detector_*', 'detection_summary.csv'};
cfg.output.own_files.decoding  = {'decoding_*'}; % <stimulus>/Decoding (stimulus decoding outputs)
cfg.output.own_files.replay    = {'replay_*'};   % <stimulus>/replay (candidate events, sequence-likelihood replay)
cfg.probe_id      = 0;      % imec probe index (single NPX probe per animal)
cfg.sorter        = 'ks4';  % loads extracted_clusters_<sorter>.mat

%% Channel selection (depths come from Ellie's depths_from_PSD / earliest_V1sink_CSD)
cfg.ca1_shank_preference = [4 3]; % most anterior CA1 shanks, first available is used
cfg.l5_shank             = [];    % [] = use the same shank as the CA1 (ripple) channel
cfg.depth_ranges.V1  = [-330 700]; % um relative to L5 depth (Ellie, after Senzai 2019)
cfg.depth_ranges.CA1 = [-150 150]; % um relative to CA1 depth
cfg.depth_ranges.L4  = [-60 60];   % um relative to L4 (CSD sink) depth

%% LFP loading
cfg.lfp.target_fs = 1250;   % Hz; detect_behavioural_states_masa assumes 1250 Hz
cfg.lfp.car       = 'auto'; % 'auto' = common average reference only if the file is not CatGT tcat (already -gblcar)
cfg.lfp.chunk_sec = 60;     % read chunk size

%% Mobility from video pixel change (Ellie's Sleep_analysis_Ellie method)
cfg.behaviour.smooth_win_sec   = 1;    % movmean window applied to Behaviour.mobility
cfg.behaviour.change_threshold = 2000; % |diff(smoothed mobility)| above this = moving
cfg.behaviour.merge_gap_sec    = 1;    % pauses shorter than this are treated as movement
% NB: FindRipples_masa / DetectSpindles_masa additionally exclude samples with
% Behaviour.mobility_zscore > 3 (hard-coded inside those detectors).

%% Sleep staging (detect_behavioural_states_masa)
cfg.sleep.speed_threshold = 1; % immobile where interpolated movement trace < 1

%% Ripples (FindRipples_masa) - onset/offset z = 1, peak z >= 5
cfg.ripple.passband    = [125 300]; % Hz
cfg.ripple.thresholds  = [1 5];     % [onset/offset, peak] z-scored ripple-band envelope
cfg.ripple.durations   = [30 200];  % [min inter-ripple interval, max duration] ms
cfg.ripple.minDuration = 30;        % ms

%% Spindles (DetectSpindles_masa, parameters as in Main_spindle_analysis_addon)
% Both detector methods are run and saved separately:
%   'new'     linear Hilbert envelope (320 ms smoothing) z-scored; + rel. power and cycle checks
%   'default' FMAToolbox FindSpindles steps: smoothed squared envelope of the sigma signal
%             z-scored against SWS (detection still runs on the whole recording), crossings
%             < 50 ms dropped, events < merge_gap apart merged, onset correction, merged
%             again; no rel. power / cycle checks (min_rel_power, min_cycles unused)
cfg.spindle.methods       = {'new', 'default'};
cfg.spindle.passband      = [9 17];    % Hz
cfg.spindle.thresholds    = struct('new', [1.5 3], 'default', [2.5 5]); % [onset/offset, peak] per method
cfg.spindle.durations     = struct('new', [400 3000], 'default', [400 3500]); % [min max] ms (FindSpindles: 400-3500)
cfg.spindle.merge_gap     = 0.4;        % s, 'default' only (Maingret et al. 2016; FindSpindles uses 0.5)
cfg.spindle.min_rel_power = 0.15;       % sigma / broadband (1-30 Hz) power ('new' only)
cfg.spindle.min_cycles    = 3;          % ('new' only)
cfg.spindle.delta_reject  = false;
cfg.spindle.delta_thresh  = 3;

%% Slow oscillations / UP-DOWN states (DetectSlowWaves_masa)
cfg.so.filterparms.deltafilter    = [0.5 4];   % Hz
cfg.so.filterparms.gammafilter    = [100 400]; % Hz
cfg.so.filterparms.gammasmoothwin = 0.08;      % s
cfg.so.filterparms.gammanormwin   = 20;        % s
cfg.so.sensitivity        = 0.6;
cfg.so.mua_quality_filter = true;  % true = only good V1 units (as Sleep_analysis_Ellie); all 400 V1 units runs out of memory in DetectSlowWaves_masa

%% Plotting
cfg.plot.visible = 'off';        % figures are saved, not shown, during batch runs
cfg.plot.svg_min_points = 5000;  % figures with more scatter/marker points are saved as .svg + raster .pdf
end
