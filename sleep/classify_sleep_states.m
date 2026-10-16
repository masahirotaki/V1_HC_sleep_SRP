function states = classify_sleep_states(ctx_lfp, hpc_lfp, tvec, Behaviour, cfg)
% CLASSIFY_SLEEP_STATES  Movement / quiet wake / SWS / REM with detect_behavioural_states_masa.
%
%   ctx_lfp  cortical (L5) LFP, 1250 Hz   -> SWS from high delta (0.5-4 Hz) power during immobility
%   hpc_lfp  CA1 LFP, 1250 Hz             -> REM from theta > delta power, preceded by SWS
%   Immobility = Behaviour.moving (video mobility) interpolated onto tvec < cfg.sleep.speed_threshold
%
% states.<SWS|REM|quietWake|movement>  [start stop] intervals in s (SpikeGLX time)

speed = interp1(Behaviour.tvec, double(Behaviour.moving), tvec, 'linear');
[~, quietWake, SWS, REM, movement] = detect_behavioural_states_masa( ...
    [tvec(:) ctx_lfp(:)], [tvec(:) hpc_lfp(:)], [tvec(:) speed(:)], cfg.sleep.speed_threshold);

states.SWS       = SWS;
states.REM       = REM;
states.quietWake = quietWake;
states.movement  = movement;
states.recording_range = tvec([1 end]);

names = {'SWS', 'REM', 'quietWake', 'movement'};
total = diff(states.recording_range);
n_epochs = cellfun(@(n) size(states.(n), 1), names)';
total_sec = cellfun(@(n) sum(diff(states.(n), [], 2)), names)';
states.summary = table(string(names'), n_epochs, total_sec, total_sec / total, ...
    'VariableNames', {'state', 'n_epochs', 'total_sec', 'fraction'});
disp(states.summary);
end
