function spindles = detect_spindles(l5_lfp, tvec, fs, Behaviour, states, channel, method, cfg)
% DETECT_SPINDLES  Cortical (L5) sleep spindles with DetectSpindles_masa.
%
% method = 'new' or 'default' (see cfg.spindle). Sigma-band (cfg.spindle.passband)
% envelope thresholded at [onset peak] = cfg.spindle.thresholds.(method):
%   'new'     filtered by duration, relative sigma power (sigma / 1-30 Hz) and
%             >= min_cycles oscillation cycles
%   'default' FMAToolbox FindSpindles steps, z-scored against the SWS epochs in
%             STATES (events are still detected in every state), events closer than
%             cfg.spindle.merge_gap merged, onset correction
% Samples with mobility_zscore > 3 are excluded.
%
% spindles.onset / .offset / .peaktimes (s), .peak_z_amp (= .peak_zscore, for the
%         shared plotting code), .peak_raw_amp (uV), .rel_power, .n_cycles,
%         .duration (s), .state, .channel, .method

p = cfg.spindle;
p.method = method;
p.thresholds = cfg.spindle.thresholds.(method);
p.durations = cfg.spindle.durations.(method);
spindles = DetectSpindles_masa(l5_lfp(:), tvec(:), 'frequency', fs, 'method', method, ...
    'passband', p.passband, 'thresholds', p.thresholds, 'durations', p.durations, ...
    'min_rel_power', p.min_rel_power, 'min_cycles', p.min_cycles, ...
    'delta_reject', p.delta_reject, 'delta_thresh', p.delta_thresh, ...
    'intervals', states.SWS, 'merge_gap', p.merge_gap, ...
    'behaviour', Behaviour, 'noise', [], 'show', 'off', 'best_channel', channel.chanmap_row);

spindles.peak_zscore = spindles.peak_z_amp;
spindles.duration = spindles.offset - spindles.onset;
spindles.state = assign_event_states(spindles.peaktimes, states);
spindles.channel = channel;
spindles.method = method;
spindles.params = rmfield(p, 'methods');
fprintf('Spindles (%s): %d detected (%d in SWS)\n', method, numel(spindles.peaktimes), ...
    sum(spindles.state == "SWS"));
end
