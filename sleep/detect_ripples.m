function ripples = detect_ripples(ca1_lfp, tvec, fs, Behaviour, states, channel, cfg)
% DETECT_RIPPLES  CA1 sharp-wave ripples with FindRipples_masa.
%
% Ripple-band (cfg.ripple.passband) Hilbert envelope is z-scored over the whole
% recording; events start/stop at z = thresholds(1), must peak above
% thresholds(2), are merged if closer than durations(1) ms, kept if between
% minDuration and durations(2) ms, and are only detected while
% Behaviour.mobility_zscore < 3.
%
% ripples.onset / .offset / .peaktimes (s; peak = trough of the filtered ripple)
%        .peak_zscore, .duration (s), .state, .channel

ripples = FindRipples_masa(ca1_lfp(:), tvec(:), 'behaviour', Behaviour, ...
    'minDuration', cfg.ripple.minDuration, 'durations', cfg.ripple.durations, ...
    'frequency', fs, 'noise', [], 'passband', cfg.ripple.passband, ...
    'thresholds', cfg.ripple.thresholds, 'show', 'off', 'best_channel', channel.chanmap_row);

ripples.duration = ripples.offset - ripples.onset;
ripples.state = assign_event_states(ripples.peaktimes, states);
ripples.channel = channel;
ripples.params = cfg.ripple;
fprintf('Ripples: %d detected (%d in SWS)\n', numel(ripples.peaktimes), sum(ripples.state == "SWS"));
end
