function so = detect_slow_oscillations(l5_lfp, tvec, units, states, channel, cfg)
% DETECT_SLOW_OSCILLATIONS  Slow oscillations and UP/DOWN states with DetectSlowWaves_masa.
%
% Within SWS, DOWN states are windows where the delta-filtered (cfg.so deltafilter,
% 0.5-4 Hz) L5 LFP has a positive peak AND high-gamma (100-400 Hz) power dips.
% Thresholds are calibrated so that V1 MUA (all V1 units, or only good ones if
% cfg.so.mua_quality_filter) is low inside DOWN states; DOWNs with too much
% spiking are rejected. UP states are the intervals between consecutive DOWNs.
%
% so.peaktimes   slow-wave (DOWN-state) peak times (s), .peak_magnitude
% so.DOWN, so.UP [start stop] intervals (s); .DOWN_duration, .UP_duration (s)
% so.detector_output   untouched DetectSlowWaves_masa output (buzcode format)

so.channel = channel;
so.params = cfg.so;

v1 = units.region == "V1";
if cfg.so.mua_quality_filter; v1 = v1 & units.is_good; end
mua.cluster_id  = units.cluster_id(v1);
in_mua = ismember(units.spike_id, mua.cluster_id);
mua.spike_times = units.spike_times(in_mua);
mua.spike_id    = units.spike_id(in_mua);
so.mua_cluster_id = mua.cluster_id;

if isempty(states.SWS)
    warning('No SWS detected: skipping slow-oscillation detection.');
    [so.peaktimes, so.peak_magnitude, so.DOWN_duration, so.UP_duration] = deal(zeros(0, 1));
    [so.DOWN, so.UP] = deal(zeros(0, 2));
    so.detector_output = [];
    return
end

args = {'filterparms', cfg.so.filterparms, 'time', tvec(:)', 'lfp', l5_lfp(:)', ...
    'NREMInts', states.SWS, 'sensitivity', cfg.so.sensitivity, 'showFig', true};
if isempty(mua.cluster_id)
    warning('No V1 units: DOWN states are not validated with spiking.');
    args = [args, {'noSpikes', true}];
else
    args = [args, {'spikes', mua}];
end
SW = DetectSlowWaves_masa(args{:});

so.peaktimes      = SW.timestamps(:);
so.peak_magnitude = SW.SWpeakmag(:);
so.DOWN           = SW.ints.DOWN;
so.UP             = SW.ints.UP;
so.DOWN_duration  = diff(so.DOWN, [], 2);
so.UP_duration    = diff(so.UP, [], 2);
so.detector_output = SW;
fprintf('Slow oscillations: %d DOWN states, %d UP states (%d V1 units for MUA)\n', ...
    size(so.DOWN, 1), size(so.UP, 1), numel(mua.cluster_id));
end
