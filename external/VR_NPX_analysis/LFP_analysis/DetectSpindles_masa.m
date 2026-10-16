function [spindles] = DetectSpindles_masa(lfp, timevec, varargin)
%DetectSpindles_masa - Detect cortical sleep spindles with improved
%   delta-wave rejection (relative sigma power + cycle count validation).
%
%   Improved spindle detection based on literature best practices:
%     - Narrower default passband (10-16 Hz) to reduce harmonic leakage
%     - Relative sigma power criterion (YASA-inspired) to reject delta harmonics
%     - Zero-crossing cycle count validation (>= 3 cycles required)
%     - Optional delta co-occurrence flagging
%     - Higher default low threshold (1.5 SD) to reduce false positives
%
% USAGE
%    [spindles] = DetectSpindles_masa(lfp, timevec, <options>)
%
% INPUTS
%    lfp            Unfiltered LFP (one channel), nsample x 1
%    timevec        Timestamps matching LFP samples
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'thresholds'      [low high] in SD (default = [1.5 3])
%     'durations'       [min max] spindle duration in ms (default = [400 3000])
%     'frequency'       Sampling rate in Hz (default = 1000)
%     'passband'        [lo hi] spindle frequency band (default = [10 16])
%     'broadband'       [lo hi] for relative power calc (default = [1 30])
%     'min_rel_power'   Minimum relative sigma power (default = 0.15)
%     'min_cycles'      Minimum oscillatory cycles (default = 3)
%     'delta_reject'    Reject events with high delta power (default = true)
%     'delta_band'      Delta frequency band (default = [0.5 4])
%     'delta_thresh'    Delta zscore threshold for rejection (default = 3)
%     'behaviour'       Behaviour struct with speed/mobility fields
%     'noise'           Noise channel LFP for artifact rejection
%     'show'            'on' or 'off' for plotting (default = 'off')
%     'best_channel'    Channel index for metadata
%     'saveMat'         Save results to .mat (default = false)
%     'savepath'        Path for saving
%     'method'          'new' or 'default' (see METHODS below)
%     'intervals'       [start stop] times (s), e.g. SWS epochs, used as the
%                       z-score reference for 'default' (default = [] = whole
%                       recording). Detection still runs on the whole recording.
%     'merge_gap'       'default' only: merge events closer than this (s)
%                       (default = 0.4, as Maingret et al. 2016; FMA uses 0.5)
%    =========================================================================
%
% METHODS
%    'new'      Linear Hilbert envelope (320 ms smoothing), z-scored; thresholds
%               in SD; relative sigma power and cycle-count criteria.
%    'default'  FMAToolbox FindSpindles steps: squared Hilbert envelope of the
%               z-scored sigma signal, 320 ms boxcar (filtfilt); crossings of
%               thresholds(1) shorter than 50 ms dropped; events closer than
%               merge_gap merged; peak >= thresholds(2); onset correction (start
%               re-thresholded at peak/3 when that exceeds thresholds(1), searching
%               back to the previous spindle's end, so a leading delta wave does
%               not pull the start earlier); merged again; durations applied.
%               No relative power or cycle criteria.
%
% OUTPUT
%    spindles       Struct with fields:
%                     .onset          Nx1 spindle start times
%                     .offset         Nx1 spindle end times
%                     .peaktimes      Nx1 peak power timestamps
%                     .peak_z_amp     Nx1 peak z-scored amplitude
%                     .peak_raw_amp   Nx1 peak raw envelope amplitude (uV)
%                     .rel_power      Nx1 mean relative sigma power
%                     .n_cycles       Nx1 number of oscillatory cycles
%                     .best_channel   Channel used for detection
%                     .detectorinfo   Struct with detection parameters
%
% REFERENCES
%   Lacourse et al. (2018) - multi-feature spindle detection (basis for YASA)
%   Latchoumane et al. (2017) PNAS - SO-spindle-ripple coupling
%   Peyrache et al. (2011) Nat Neurosci - Hilbert-based spindle detection
%   FMAToolbox FindSpindles (Zugaro, Maingret, Todorova) - 'default' method
%   Maingret et al. (2016) Nat Neurosci - 0.4 s merging
%
% Modified from FindSpindles_masa.m by Masahiro Takigawa
% Masahiro Takigawa, 2026

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
% modified by Masahiro Takigawa, 2023, 2026

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% Parse inputs
p = inputParser;
addParameter(p, 'thresholds',    [1.5 3],     @isnumeric);
addParameter(p, 'durations',     [400 3000],   @isnumeric);
addParameter(p, 'frequency',     1000,         @isnumeric);
addParameter(p, 'passband',      [9 17],      @isnumeric);
addParameter(p, 'broadband',     [1 30],       @isnumeric);
addParameter(p, 'min_rel_power', 0.15,         @isnumeric);
addParameter(p, 'min_cycles',    3,            @isnumeric);
addParameter(p, 'delta_reject',  false,         @islogical);
addParameter(p, 'delta_band',    [0.5 4],      @isnumeric);
addParameter(p, 'delta_thresh',  3,            @isnumeric);
addParameter(p, 'behaviour',     [],           @isstruct);
addParameter(p, 'noise',         [],           @ismatrix);
addParameter(p, 'show',          'off',        @isstr);
addParameter(p, 'best_channel',  [],           @isnumeric);
addParameter(p, 'saveMat',       false,        @islogical);
addParameter(p, 'savepath',      [],           @isstr);
addParameter(p, 'method',      [],           @isstr);
addParameter(p, 'intervals',     [],           @isnumeric);
addParameter(p, 'merge_gap',     0.4,          @isnumeric);
parse(p, varargin{:});

frequency       = p.Results.frequency;
passband        = p.Results.passband;
broadband_range = p.Results.broadband;
lowThresh       = p.Results.thresholds(1);
highThresh      = p.Results.thresholds(2);
minDuration     = p.Results.durations(1) / 1000;  % convert ms -> s
maxDuration     = p.Results.durations(2) / 1000;
min_rel_power   = p.Results.min_rel_power;
min_cycles      = p.Results.min_cycles;
delta_reject    = p.Results.delta_reject;
delta_band      = p.Results.delta_band;
delta_thresh    = p.Results.delta_thresh;
behaviour       = p.Results.behaviour;
noise           = p.Results.noise;
show            = p.Results.show;
best_channel    = p.Results.best_channel;
method          = p.Results.method;
intervals       = p.Results.intervals;
merge_gap       = p.Results.merge_gap;

if strcmp(method, 'default') && ismember('thresholds', p.UsingDefaults)
    lowThresh = 2.5;
    highThresh = 5;
end
if strcmp(method, 'default') && ismember('durations', p.UsingDefaults)
    maxDuration = 3.5;   % FindSpindles default [400 3500] ms
end

% Ensure column vectors
lfp = lfp(:);
timevec = timevec(:);
lfp = lfp-mean(lfp);


%% 1. Bandpass filter in sigma band
filter_order = round(6 * frequency / (max(passband) - min(passband)));
norm_freq = passband / (frequency / 2);
b_sigma = fir1(filter_order, norm_freq, 'bandpass');
sigma_signal = filtfilt(b_sigma, 1, lfp);

%% 2. Broadband filter (1-30 Hz) for relative power calculation
% filter_order_bb = round(6 * frequency / (max(broadband_range) - min(broadband_range)));
% norm_freq_bb = broadband_range / (frequency / 2);
% b_broadband = fir1(filter_order_bb, norm_freq_bb, 'bandpass');
% broadband_signal = filtfilt(b_broadband, 1, lfp);

% Design a 3rd-order Butterworth bandpass filter
% (Note: filtfilt doubles the order to 6, which is ideal for EEG)
norm_freq_bb = broadband_range / (frequency / 2);
[b_broadband, a_broadband] = butter(3, norm_freq_bb, 'bandpass');

% Apply the filter to the raw lfp 
% (The sharp 1 Hz cutoff will automatically kill the DC offset)
broadband_signal = filtfilt(b_broadband, a_broadband, lfp);

%% 3. Compute envelopes
% rms_window = round(frequency / 5);  % ~200 ms window

% Z-score the filtered signal
% 1. Raw Hilbert Envelopes for Relative Power
% Get raw squared amplitude (microvolts squared)
raw_sq_sigma = abs(hilbert(sigma_signal)).^2;
raw_sq_broad = abs(hilbert(broadband_signal)).^2;

% Smooth both exactly the same way (e.g. 100ms or 320ms window)
if strcmp(method, 'new')
    % 100 ms moving average (1 cycle) is optimal for [1.5, 3] thresholds in mice
    windowLength = round(frequency * 0.320);
else
    % Classic FMA FindSpindles uses exactly 320 ms
    windowLength = round(frequency * 0.320);
end
smooth_win = ones(windowLength, 1) / windowLength;

power_sigma = filtfilt(smooth_win, 1, raw_sq_sigma);
power_broad = filtfilt(smooth_win, 1, raw_sq_broad);

% Compute relative sigma power (YASA-inspired criterion)
power_broad(power_broad < eps) = eps;
relative_sigma_power = power_sigma ./ power_broad;

if strcmp(method, 'new')
    % --- NEW HYBRID METHOD ---

    %Calculate Linear Z-scored Envelope for [1.5, 3] Thresholding
    % We take the linear amplitude, smooth it, THEN z-score it
    % This preserves the meaning of your [1.5, 3] SD thresholds!
    sigma_envelope = abs(hilbert(sigma_signal));
    sigma_envelope = filtfilt(smooth_win, 1, sigma_envelope);
    zscored_sigma  = zscore(sigma_envelope);

    % (Optional) Keep broadband linear envelope for downstream plotting compatibility
    broadband_envelope = sqrt(power_broad);

else % 'default' (classic FMA FindSpindles logic)

    % 2. Normalized Hilbert Envelope strictly for Thresholding
    % Z-score the raw filtered signal first (Classic Buzsaki method). FindSpindles
    % expects SWS-only input, so the mean/SD can come from 'intervals' (e.g. SWS).
    ref = true(size(sigma_signal));
    if ~isempty(intervals)
        ref = false(size(sigma_signal));
        for k = 1:size(intervals, 1)
            ref(timevec >= intervals(k, 1) & timevec <= intervals(k, 2)) = true;
        end
    end
    z_sigma = (sigma_signal - mean(sigma_signal(ref))) / std(sigma_signal(ref));

    % Get squared envelope of the Z-scored signal
    z_sq_sigma = abs(hilbert(z_sigma)).^2;

    % Smooth it
    zscored_sigma = filtfilt(smooth_win, 1, z_sq_sigma);

    % Now you can apply the strict [2.5, 5] Buzsaki thresholds to `zscored_sigma_envelope`!
end




%% 5. Optional: Delta power for co-occurrence check
if delta_reject
    filter_order_delta = round(6 * frequency / (max(delta_band) - min(delta_band)));
    norm_freq_delta = delta_band / (frequency / 2);
    b_delta = fir1(filter_order_delta, norm_freq_delta, 'bandpass');
    delta_signal = filtfilt(b_delta, 1, lfp);
    delta_envelope = abs(hilbert(delta_signal));
    zscored_delta = zscore(delta_envelope);
end

%% 6. Speed exclusion mask
speed_mask = true(size(timevec));
if ~isempty(behaviour)
    if isfield(behaviour, 'mobility_zscore')
        speed_threshold = 3;  % zscore
        speed = interp1(behaviour.sglxTime, behaviour.mobility_zscore, timevec, 'nearest');
    elseif isfield(behaviour, 'speed')
        speed_threshold = 5;  % cm/s
        speed = interp1(behaviour.sglxTime, behaviour.speed, timevec, 'nearest');
    else
        speed = zeros(size(timevec));
        speed_threshold = inf;
    end
    speed(isnan(speed)) = 0;
    speed_mask = speed < speed_threshold;
end

%% 7. Initial threshold detection
thresholded = (zscored_sigma > lowThresh) & speed_mask;

% Find start/stop indices
d = diff([0; thresholded; 0]);
start_idx = find(d > 0);
stop_idx  = find(d < 0) - 1;

if isempty(start_idx) || isempty(stop_idx)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindle candidates found after initial thresholding.');
    return
end

% Ensure equal length
n_events = min(length(start_idx), length(stop_idx));
start_idx = start_idx(1:n_events);
stop_idx = stop_idx(1:n_events);

firstPass = [start_idx, stop_idx];
disp(['DetectSpindles: ' num2str(size(firstPass, 1)) ' candidates after thresholding.']);

if strcmp(method, 'new')
    %% 8. Peak threshold filter
    secondPass = [];
    peakNormedPower = [];
    for i = 1:size(firstPass, 1)
        segment = zscored_sigma(firstPass(i,1):firstPass(i,2));
        [maxVal, ~] = max(segment);
        if maxVal > highThresh
            secondPass = [secondPass; firstPass(i,:)];
            peakNormedPower = [peakNormedPower; maxVal];
        end
    end

    if isempty(secondPass)
        spindles = make_empty_spindles(best_channel, p.Results);
        disp('No spindles survived peak threshold.');
        return
    end
    disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' candidates after peak threshold.']);

    %% 9. Duration filter
    durations_s = timevec(secondPass(:,2)) - timevec(secondPass(:,1));
    valid_dur = (durations_s >= minDuration) & (durations_s <= maxDuration);
    secondPass = secondPass(valid_dur, :);
    peakNormedPower = peakNormedPower(valid_dur);
else
    %% 8-9. FMAToolbox FindSpindles steps: drop < 50 ms bouts, merge, peak test,
    % onset correction, merge again, duration filter
    [secondPass, peakNormedPower] = findspindles_steps(zscored_sigma, firstPass, ...
        frequency, lowThresh, highThresh, merge_gap, minDuration, maxDuration);
end

if isempty(secondPass)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindles survived duration filter.');
    return
end
disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' candidates after duration filter.']);

%% 10. Relative sigma power filter (key improvement over FindSpindles_masa)
mean_rel_power = zeros(size(secondPass, 1), 1);
for i = 1:size(secondPass, 1)
    mean_rel_power(i) = mean(relative_sigma_power(secondPass(i,1):secondPass(i,2)));
end

if strcmp(method, 'new')
    valid_rel = mean_rel_power >= min_rel_power;
else
    valid_rel = true(size(secondPass, 1), 1); % Skip relative power check in default mode
end
secondPass = secondPass(valid_rel, :);
peakNormedPower = peakNormedPower(valid_rel);
mean_rel_power = mean_rel_power(valid_rel);

if isempty(secondPass)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindles survived relative power filter.');
    return
end
disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' candidates after relative sigma power filter.']);

%% 11. Zero-crossing cycle count validation
n_cycles_vec = zeros(size(secondPass, 1), 1);
for i = 1:size(secondPass, 1)
    seg = sigma_signal(secondPass(i,1):secondPass(i,2));
    % Count zero crossings: each pair = 1 cycle
    zc = sum(diff(sign(seg)) ~= 0);
    n_cycles_vec(i) = zc / 2;  % half-cycles -> full cycles
end

if strcmp(method, 'new')
    valid_cycles = n_cycles_vec >= min_cycles;
else
    valid_cycles = true(size(secondPass, 1), 1); % Skip cycle check in default mode
end
secondPass = secondPass(valid_cycles, :);
peakNormedPower = peakNormedPower(valid_cycles);
mean_rel_power = mean_rel_power(valid_cycles);
n_cycles_vec = n_cycles_vec(valid_cycles);

if isempty(secondPass)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindles survived cycle count filter.');
    return
end
disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' candidates after cycle count filter.']);

%% 12. Delta co-occurrence rejection
delta_rejected = false(size(secondPass, 1), 1);
if delta_reject
    for i = 1:size(secondPass, 1)
        mean_delta = mean(zscored_delta(secondPass(i,1):secondPass(i,2)));
        if mean_delta > delta_thresh
            delta_rejected(i) = true;
        end
    end
    
    n_delta_rejected = sum(delta_rejected);
    secondPass = secondPass(~delta_rejected, :);
    peakNormedPower = peakNormedPower(~delta_rejected);
    mean_rel_power = mean_rel_power(~delta_rejected);
    n_cycles_vec = n_cycles_vec(~delta_rejected);
    
    if n_delta_rejected > 0
        disp(['DetectSpindles: Rejected ' num2str(n_delta_rejected) ' events due to high delta co-occurrence.']);
    end
end

if isempty(secondPass)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindles survived delta rejection.');
    return
end

%% 13. Noise channel rejection
if ~isempty(noise)
    noise_filtered = filtfilt(b_sigma, 1, noise(:));
    zscored_noise = zscore(abs(hilbert(noise_filtered)));
    
    excluded = false(size(secondPass, 1), 1);
    for i = 1:size(secondPass, 1)
        if any(zscored_noise(secondPass(i,1):secondPass(i,2)) > 5)
            excluded(i) = true;
        end
    end
    
    secondPass = secondPass(~excluded, :);
    peakNormedPower = peakNormedPower(~excluded);
    mean_rel_power = mean_rel_power(~excluded);
    n_cycles_vec = n_cycles_vec(~excluded);
    disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' candidates after noise channel rejection.']);
end

if isempty(secondPass)
    spindles = make_empty_spindles(best_channel, p.Results);
    disp('No spindles survived noise rejection.');
    return
end

%% 14. Find peak positions (negative peak in filtered signal)
peakPosition = zeros(size(secondPass, 1), 1);
for i = 1:size(secondPass, 1)
    [~, minIdx] = min(sigma_signal(secondPass(i,1):secondPass(i,2)));
    peakPosition(i) = minIdx + secondPass(i,1) - 1;
end

%% 15. Recalculate peak power and extract aligned amplitudes
% Compute the full raw instantaneous envelope once
raw_sigma_envelope = abs(hilbert(sigma_signal)); 

raw_hilbert_amplitude = zeros(size(secondPass, 1), 1);

for i = 1:size(secondPass, 1)
    seg_start = secondPass(i,1);
    seg_end   = secondPass(i,2);
    
    % 1. Find the exact index where the Z-scored detection envelope peaks
    [maxVal, maxIdx] = max(zscored_sigma(seg_start:seg_end));
    
    % 2. Save the Z-scored amplitude (this overwrites peakNormedPower)
    peakNormedPower(i) = maxVal;
    
    % 3. Save the Raw amplitude at that EXACT same index
    raw_hilbert_amplitude(i) = raw_sigma_envelope(seg_start + maxIdx - 1);
end

disp(['DetectSpindles: ' num2str(size(secondPass, 1)) ' spindle events detected (final).']);

%% 16. Build output struct
spindles.onset         = timevec(secondPass(:, 1));
spindles.offset        = timevec(secondPass(:, 2));
spindles.peaktimes     = timevec(peakPosition);
spindles.peak_z_amp    = peakNormedPower;         
spindles.peak_raw_amp  = raw_hilbert_amplitude;   
spindles.rel_power     = mean_rel_power;
spindles.n_cycles      = n_cycles_vec;
spindles.best_channel  = best_channel;

% Detector info
detectorinfo.detectorname    = 'DetectSpindles_masa';
detectorinfo.detectiondate   = datetime('today');
detectorinfo.detectionparms  = p.Results;
if isfield(detectorinfo.detectionparms, 'noise')
    detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms, 'noise');
end
if isfield(detectorinfo.detectionparms, 'behaviour')
    detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms, 'behaviour');
end
spindles.detectorinfo = detectorinfo;

%% 17. Plotting
if strcmp(show, 'on')
    plot_spindle_detection(timevec, lfp, sigma_signal, zscored_sigma, ...
        relative_sigma_power, spindles, lowThresh, highThresh, ...
        min_rel_power, behaviour, frequency, passband, method);
end

%% 18. Save
% if p.Results.saveMat
%     save(fullfile(p.Results.savepath, 'detected_spindle_events.mat'), 'spindles');
% end

end


%% ========================================================================
%  Helper functions
%  ========================================================================

function [events, peaks] = findspindles_steps(z, firstPass, fs, thr, minPeak, merge_gap, minDur, maxDur)
%FINDSPINDLES_STEPS Event logic of FMAToolbox FindSpindles on a thresholded envelope.
%   z          smoothed squared envelope of the z-scored sigma signal
%   firstPass  [start stop] sample indices of crossings of thr
    on = firstPass(:, 1); off = firstPass(:, 2);
    keep = (off - on + 1) / fs >= 0.050;          % drop bouts < 50 ms
    on = on(keep); off = off(keep);
    gap = round(merge_gap * fs);
    [on, off] = merge_close(on, off, zeros(size(on)), zeros(size(on)), gap);
    pv = zeros(size(on)); pk = zeros(size(on));
    for i = 1:numel(on)
        [pv(i), im] = max(z(on(i):off(i)));
        pk(i) = on(i) + im - 1;
    end
    keep = pv >= minPeak;
    on = on(keep); off = off(keep); pv = pv(keep); pk = pk(keep);
    % onset correction: a leading delta wave can pull the start earlier, so for
    % strong events the start is re-thresholded at peak/3
    for i = 1:numel(on)
        newThr = pv(i) / 3;
        if newThr > thr
            if i == 1; lo = 1; else; lo = off(i - 1) + 1; end
            j = find(z(lo:pk(i)) < newThr, 1, 'last');
            if ~isempty(j); on(i) = lo + j; end   % first sample above newThr
        end
    end
    [on, off, pv] = merge_close(on, off, pv, pk, gap);
    dur = (off - on) / fs;
    keep = dur >= minDur & dur <= maxDur;
    events = [on(keep) off(keep)];
    peaks = pv(keep);
end


function [on, off, pv, pk] = merge_close(on, off, pv, pk, gap)
%MERGE_CLOSE Merge events separated by fewer than GAP samples (larger peak kept).
    i = 1;
    while i < numel(on)
        if on(i + 1) - off(i) < gap
            off(i) = off(i + 1);
            if pv(i + 1) > pv(i); pv(i) = pv(i + 1); pk(i) = pk(i + 1); end
            on(i + 1) = []; off(i + 1) = []; pv(i + 1) = []; pk(i + 1) = [];
        else
            i = i + 1;
        end
    end
end


function spindles = make_empty_spindles(best_channel, params)
%MAKE_EMPTY_SPINDLES Create an empty spindles struct with consistent fields.
    spindles.onset       = [];
    spindles.offset      = [];
    spindles.peaktimes   = [];
    spindles.peak_z_amp  = [];
    spindles.peak_raw_amp = [];
    spindles.rel_power   = [];
    spindles.n_cycles    = [];
    spindles.best_channel = best_channel;
    
    detectorinfo.detectorname   = 'DetectSpindles_masa';
    detectorinfo.detectiondate  = datetime('today');
    detectorinfo.detectionparms = params;
    if isfield(detectorinfo.detectionparms, 'noise')
        detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms, 'noise');
    end
    if isfield(detectorinfo.detectionparms, 'behaviour')
        detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms, 'behaviour');
    end
    spindles.detectorinfo = detectorinfo;
end


function plot_spindle_detection(timevec, lfp, sigma_signal, zscored_sigma, ...
    relative_sigma_power, spindles, lowThresh, highThresh, ...
    min_rel_power, behaviour, frequency, passband, method)
%PLOT_SPINDLE_DETECTION Visualize detected spindles.

    if isempty(spindles.onset) || length(spindles.onset) < 3
        disp('Too few spindles to plot.');
        return
    end
    
    % Pick a segment with spindles in the middle of the recording
    mid_idx = round(length(spindles.onset) / 2);
    first_sp = max(1, mid_idx - 2);
    last_sp  = min(length(spindles.onset), mid_idx + 3);
    
    t_start = spindles.onset(first_sp) - 2;
    t_end   = spindles.offset(last_sp) + 2;
    time_mask = timevec >= t_start & timevec <= t_end;
    
    figure('Name', sprintf('spindle_detection_%s', method), ...
           'Position', [100 100 1200 800]);
    
    % Panel 1: Raw LFP with spindle-band overlay
    ax1 = subplot(4, 1, 1);
    plot(timevec(time_mask), lfp(time_mask), 'Color', [0.5 0.5 0.5]);
    hold on;
    plot(timevec(time_mask), sigma_signal(time_mask), 'b', 'LineWidth', 1);
    for j = first_sp:last_sp
        xline(spindles.onset(j), 'g', 'LineWidth', 1.5);
        xline(spindles.offset(j), 'r', 'LineWidth', 1.5);
    end
    title(sprintf('Raw LFP (grey) + Sigma filtered (%d-%d Hz, blue)', passband(1), passband(2)));
    ylabel('Amplitude');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    % Panel 2: Z-scored sigma power with thresholds
    ax2 = subplot(4, 1, 2);
    plot(timevec(time_mask), zscored_sigma(time_mask), 'k');
    hold on;
    yline(lowThresh, '--', 'Low threshold', 'Color', [0.2 0.6 0.2]);
    yline(highThresh, '--', 'High threshold', 'Color', [0.8 0.2 0.2]);
    for j = first_sp:last_sp
        xline(spindles.onset(j), 'g', 'LineWidth', 1.5);
        xline(spindles.offset(j), 'r', 'LineWidth', 1.5);
    end
    title('Z-scored sigma RMS power');
    ylabel('Z-score');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    % Panel 3: Relative sigma power
    ax3 = subplot(4, 1, 3);
    plot(timevec(time_mask), relative_sigma_power(time_mask), 'Color', [0.8 0.4 0]);
    hold on;
    yline(min_rel_power, '--', 'Min relative power', 'Color', [0.6 0.2 0]);
    for j = first_sp:last_sp
        xline(spindles.onset(j), 'g', 'LineWidth', 1.5);
        xline(spindles.offset(j), 'r', 'LineWidth', 1.5);
    end
    title('Relative sigma power (sigma^2 / broadband^2)');
    ylabel('Ratio');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    % Panel 4: Speed (if available)
    ax4 = subplot(4, 1, 4);
    if ~isempty(behaviour)
        if isfield(behaviour, 'mobility_zscore')
            spd = interp1(behaviour.sglxTime, behaviour.mobility_zscore, timevec, 'nearest');
            ylabel('Mobility (z-score)');
        elseif isfield(behaviour, 'speed')
            spd = interp1(behaviour.sglxTime, behaviour.speed, timevec, 'nearest');
            ylabel('Speed (cm/s)');
        else
            spd = zeros(size(timevec));
            ylabel('Speed');
        end
        plot(timevec(time_mask), spd(time_mask), 'Color', [0.3 0.3 0.8]);
    end
    title('Movement');
    xlabel('Time (s)');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    linkaxes([ax1 ax2 ax3 ax4], 'x');
    
    % Summary figure: spindle statistics
    figure('Name', sprintf('spindle_detection_summary_%s', method), 'Position', [100 100 1000 400]);
    
    subplot(1, 3, 1);
    histogram(spindles.peak_z_amp, 30, 'FaceColor', [0.2 0.5 0.8]);
    xlabel('Peak z-score'); ylabel('Count');
    title('Peak sigma power'); 
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    subplot(1, 3, 2);
    durations_ms = (spindles.offset - spindles.onset) * 1000;
    histogram(durations_ms, 30, 'FaceColor', [0.8 0.4 0.2]);
    xlabel('Duration (ms)'); ylabel('Count');
    title('Spindle duration');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
    
    subplot(1, 3, 3);
    histogram(spindles.rel_power, 30, 'FaceColor', [0.2 0.7 0.3]);
    xlabel('Relative sigma power'); ylabel('Count');
    title('Relative power');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 11);
end
