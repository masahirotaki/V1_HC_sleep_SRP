function results = run_sleep_detection(subject, date, stimulus_name, cfg)
% RUN_SLEEP_DETECTION  Sleep staging + ripple, spindle and slow-oscillation (UP/DOWN)
% detection for one sleep recording.
%
%   results = run_sleep_detection('M00069', '20250929', 'Sleep_Box_1', project_config())
%
% Writes into the recording's analysis folder, SUBJECTS/<subject>/analysis/<date>/<stimulus>/
% (replacing Ellie's earlier behavioural_state.mat; no other existing file is touched):
%   behavioural_state.mat        behavioural_state (SWS/REM/quietWake/movement intervals)
%   extracted_LFP.mat            LFP: CA1 / L5 / L4 channel of every shank, 1250 Hz, uV; session
%   behavioural_state.fig/.pdf   overview: mobility, spectrograms, delta, theta/delta
%   sleep_stages.fig/.pdf        hypnogram (Movement / Quiet Wake / SWS / REM rows)
%   detection/extracted_ripple_events.mat            ripples
%   detection/extracted_spindle_events_<method>.mat  spindles, one file per cfg.spindle.methods
%   detection/extracted_slow_oscillation_events.mat  slow_oscillations (UP/DOWN)
%   detection/*.fig/.pdf(/.svg), detection/detection_summary.csv

run_timer = tic;
default_visible = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', cfg.plot.visible);
restore_visible = onCleanup(@() set(groot, 'DefaultFigureVisible', default_visible));

paths = get_session_paths(subject, date, stimulus_name, cfg);
% check_output_path(paths.detection_dir, cfg);
title_str = sprintf('%s %s %s', subject, paths.date, stimulus_name);
fprintf('\n===== %s =====\n', title_str);

%% Load
layout = get_probe_layout(paths, cfg);
Behaviour = load_session_behaviour(paths, cfg);
units = load_session_clusters(paths, layout, cfg);

rows = [[layout.shank.CA1_channel]; [layout.shank.L5_channel]; [layout.shank.L4_channel]];
labels = ["CA1"; "L5"; "L4"] + "_shank" + string([layout.shank.id]);
valid = ~isnan(rows(:));
LFP = load_npx_lfp(paths.lf_bin_file, layout.chanMap(rows(valid)), cfg);
LFP.chanmap_rows = rows(valid)';
LFP.labels = labels(valid)';
LFP.shank = layout.kcoords(LFP.chanmap_rows)';
LFP.depth = layout.ycoords(LFP.chanmap_rows)';

ca1_channel = channel_info(layout, layout.ripple_channel);
l5_channel  = channel_info(layout, layout.l5_channel);
ca1 = LFP.data(find(LFP.chanmap_rows == layout.ripple_channel, 1), :);
l5  = LFP.data(find(LFP.chanmap_rows == layout.l5_channel, 1), :);
tvec = LFP.tvec;

%% Detect
states = classify_sleep_states(l5, ca1, tvec, Behaviour, cfg);
ripples = detect_ripples(ca1, tvec, LFP.fs, Behaviour, states, ca1_channel, cfg);
spindles = struct(); % one field per DetectSpindles_masa method ('new', 'default')
for m = string(cfg.spindle.methods)
    spindles.(m) = detect_spindles(l5, tvec, LFP.fs, Behaviour, states, l5_channel, char(m), cfg);
end
fig_numbers = @(figs) arrayfun(@(f) f.Number, figs); % safe when no figure is open
figs_before = fig_numbers(findall(groot, 'Type', 'figure'));
slow_oscillations = detect_slow_oscillations(l5, tvec, units, states, l5_channel, cfg);
detector_figs = findall(groot, 'Type', 'figure');
detector_figs = detector_figs(~ismember(fig_numbers(detector_figs), figs_before));
[~, order] = sort(fig_numbers(detector_figs));
detector_figs = detector_figs(order);

%% Save
% Mobility is not saved: it is recomputed from Ellie's extracted_behaviour.mat by
% load_session_behaviour. The channels / thresholds behind the states are kept in .info.
behavioural_state = states;
behavioural_state.info = struct('cortical_channel', l5_channel, 'hippocampal_channel', ca1_channel, ...
    'speed_threshold', cfg.sleep.speed_threshold, 'mobility', cfg.behaviour, 'created', datetime('now'));
backup_ellie_behavioural_state(paths.output_dir, cfg);
safe_save(fullfile(paths.output_dir, 'behavioural_state.mat'), ...
    struct('behavioural_state', behavioural_state), cfg);

session = struct('subject', subject, 'date', paths.date, 'stimulus', stimulus_name, ...
    'paths', paths, 'layout', layout, 'cfg', cfg, 'created', datetime('now'));
LFP.data = single(LFP.data);
safe_save(fullfile(paths.output_dir, 'extracted_LFP.mat'), struct('LFP', LFP, 'session', session), cfg);

% each event struct carries its own .channel and .params
safe_save(fullfile(paths.detection_dir, 'extracted_ripple_events.mat'), struct('ripples', ripples), cfg);
for m = string(fieldnames(spindles))'
    safe_save(fullfile(paths.detection_dir, sprintf('extracted_spindle_events_%s.mat', m)), ...
        struct('spindles', spindles.(m)), cfg);
end
safe_save(fullfile(paths.detection_dir, 'extracted_slow_oscillation_events.mat'), ...
    struct('slow_oscillations', slow_oscillations), cfg);

summary = detection_summary(states, ripples, spindles, slow_oscillations);
summary_file = fullfile(paths.detection_dir, 'detection_summary.csv');
check_output_path(summary_file, cfg);
check_overwrite(summary_file, cfg);
writetable(summary, summary_file);
disp(summary);

%% Plot (save_all_figures: .fig + .pdf, or .svg + raster .pdf for dense scatter plots)
% DetectSlowWaves_masa threshold-calibration figures first, so each save sees only its own figure
save_figures(detector_figs, compose('slow_oscillation_detector_%d', 1:numel(detector_figs)), ...
    paths.detection_dir, cfg);

fig = plot_sleep_states(states, Behaviour, l5, ca1, tvec, LFP.fs, cfg, title_str);
save_figures(fig, {'behavioural_state'}, paths.output_dir, cfg);
fig = plot_sleep_hypnogram(states, cfg, title_str);
save_figures(fig, {'sleep_stages'}, paths.output_dir, cfg);

opt = struct('name', 'Ripple', 'band', cfg.ripple.passband, 'avg_win', [-0.15 0.15], ...
    'example_win', [-0.1 0.1], 'peak_threshold', cfg.ripple.thresholds(2));
fig = plot_oscillation_events(ripples, ca1, tvec, LFP.fs, states, opt, cfg, title_str);
save_figures(fig, {'extracted_ripple_events'}, paths.detection_dir, cfg);

for m = string(fieldnames(spindles))'
    opt = struct('name', 'Spindle', 'band', cfg.spindle.passband, 'avg_win', [-1.5 1.5], ...
        'example_win', [-1.5 1.5], 'peak_threshold', cfg.spindle.thresholds.(m)(2));
    fig = plot_oscillation_events(spindles.(m), l5, tvec, LFP.fs, states, opt, cfg, ...
        sprintf('%s (%s method)', title_str, m));
    save_figures(fig, {sprintf('extracted_spindle_events_%s', m)}, paths.detection_dir, cfg);
end

fig = plot_slow_oscillation_events(slow_oscillations, l5, tvec, LFP.fs, units, states, cfg, title_str);
save_figures(fig, {'extracted_slow_oscillation_events'}, paths.detection_dir, cfg);

results = struct('paths', paths, 'layout', layout, 'states', states, 'ripples', ripples, ...
    'spindles', spindles, 'slow_oscillations', slow_oscillations, 'summary', summary);
fprintf('Done %s in %.1f min. Outputs in %s\n', title_str, toc(run_timer) / 60, paths.output_dir);
end


function ch = channel_info(layout, row)
ch.chanmap_row = row;
ch.bin_channel = layout.chanMap(row);
ch.shank = layout.kcoords(row);
ch.depth = layout.ycoords(row);
ch.xcoord = layout.xcoords(row);
end


function T = detection_summary(states, ripples, spindles, so)
names = ["SWS"; "REM"; "quietWake"; "movement"];
minutes = arrayfun(@(n) sum(diff(states.(n), [], 2)) / 60, names);
n_ripples  = arrayfun(@(n) sum(ripples.state == n), names);
T = table(names, minutes, n_ripples, n_ripples ./ minutes, ...
    'VariableNames', {'state', 'minutes', 'n_ripples', 'ripples_per_min'});
for m = string(fieldnames(spindles))' % one column pair per spindle method
    n_spindles = arrayfun(@(n) sum(spindles.(m).state == n), names);
    T.("n_spindles_" + m) = n_spindles;
    T.("spindles_per_min_" + m) = n_spindles ./ minutes;
end
n_so = [numel(so.peaktimes); 0; 0; 0]; % slow oscillations are detected within SWS only
T.n_slow_oscillations = n_so;
T.slow_oscillations_per_min = n_so ./ minutes;
end
