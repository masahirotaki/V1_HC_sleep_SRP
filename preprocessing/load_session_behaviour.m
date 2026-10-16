function Behaviour = load_session_behaviour(paths, cfg)
% LOAD_SESSION_BEHAVIOUR  Load Ellie's extracted_behaviour.mat and derive mobility.
%
% Behaviour (from import_and_align_Bonsai_Sleep_Ellie) is resampled to 60 Hz on
% SpikeGLX time:
%   .sglxTime / .tvec  time (s), same clock as LFP and spike times
%   .mobility          video pixel-change from Bonsai (arbitrary units)
%   .X, .Y             tracked position (pixels)
% Added here (Ellie's Sleep_analysis_Ellie method):
%   .mobility_change   |diff(movmean(mobility, smooth_win_sec))|
%   .moving            logical, mobility_change > change_threshold, with pauses
%                      shorter than merge_gap_sec filled in as movement
%   .mobility_zscore   z-scored mobility_change (FindRipples_masa / DetectSpindles_masa
%                      exclude samples with mobility_zscore > 3)

S = load(fullfile(paths.analysis_dir, 'extracted_behaviour.mat'), 'Behaviour');
Behaviour = S.Behaviour;

t = Behaviour.tvec;
fs = 1 / median(diff(t));
mobility = fillmissing(double(Behaviour.mobility), 'linear');
smoothed = movmean(mobility, max(1, round(cfg.behaviour.smooth_win_sec * fs)));
change = abs([0 diff(smoothed)]);

moving = change > cfg.behaviour.change_threshold;
idx = find(moving);
for k = find(diff(t(idx)) < cfg.behaviour.merge_gap_sec)
    moving(idx(k):idx(k + 1)) = true;
end

Behaviour.mobility_change = change;
Behaviour.moving = moving;
Behaviour.mobility_zscore = zscore(change);

fprintf('Behaviour: %.1f min, moving %.1f%% of the time\n', (t(end) - t(1)) / 60, 100 * mean(moving));
end
