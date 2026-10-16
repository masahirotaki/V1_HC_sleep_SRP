function [ripples] = FindRipples_masa(lfp,timevec,varargin)
%FindRipples - Find hippocampal ripples (100~200Hz oscillations).
%
% USAGE
%    [ripples] = bz_FindRipples(lfp.data,lfp.timestamps,<options>)
%    OR
%    [ripples] = bz_FindRipples(basepath,channel,<options>)
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria
%
% INPUTS
%    lfp            unfiltered LFP (one channel) to use (nsample x 1)
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see table below)
%
%    OR
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100]).
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization (default = all)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection
%                   (default = [130 200])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'plotType'   1=original version (several plots); 2=only raw lfp
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector
% SEE ALSO
%
%    See also bz_Filter, bz_RippleStats, bz_SaveRippleEvents, bz_PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
% modified by Masahiro Takigawa, 2023

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Based on Tingley and Buzsáki 2020
warning('this function is under development and may not work... yet')

% Default values
p = inputParser;
addParameter(p,'thresholds',[3 5],@isnumeric)
addParameter(p,'durations',[30 100],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'frequency',1000,@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'show','on',@isstr)
addParameter(p,'noise',[],@ismatrix)
addParameter(p,'passband',[125 300],@isnumeric)
addParameter(p,'speed_threshold',5,@isnumeric);
addParameter(p,'behaviour',[],@isstruct)
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'minDuration',20,@isnumeric)
addParameter(p,'plotType',1,@isnumeric)
addParameter(p,'savepath',[],@isstr)
addParameter(p,'best_channel',[],@isnumeric)
% lfp;
% timevec;

% assign parameters (either defaults or given)
parse(p,varargin{:});
frequency = p.Results.frequency;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);
minInterRippleInterval = p.Results.durations(1);
maxRippleDuration = p.Results.durations(2);
minRippleDuration = p.Results.minDuration;
plotType = p.Results.plotType;
passband = p.Results.passband;
speed_threshold = p.Results.speed_threshold;
behaviour = p.Results.behaviour;
savepath = p.Results.savepath;
best_channel = p.Results.best_channel;


filter_type  = 'bandpass';
filter_order = round(6*frequency/(max(passband)-min(passband)));  % creates filter for ripple
norm_freq_range = passband/(frequency/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
b_ripple = fir1(filter_order, norm_freq_range,filter_type);
signal = filtfilt(b_ripple,1,lfp);
zscored_ripple = zscore(abs(hilbert(signal)));

% signal = bz_Filter(double(lfp),'filter','butter','passband',passband,'order', 3);


%% Detect ripple and calculate noise

%
% % Parameters
% windowLength = frequency/frequency*11;
%
% % Square and normalize signal
% squaredSignal = signal.^2;
% % squaredSignal = abs(opsignal);
% window = ones(windowLength,1)/windowLength;


% keep = [];
% if ~isempty(restrict)
%     for i=1:size(restrict,1)
%         keep = InIntervals(timevec,restrict);
%     end
% end
% keep = logical(keep);
% zscored_ripple = zscored_ripple(keep);
% timevec = timevec(keep);
%


% [zscored_ripple,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);


% Detect ripple periods by thresholding zscored ripple power
if isfield(behaviour,'mobility_zscore')
    speed_threshold = 3; % zscore of three in terms of movement
    behaviour.speed=behaviour.mobility_zscore;
else
    speed_threshold = 5; %cm
end

speed = interp1(behaviour.sglxTime,behaviour.speed,timevec,'nearest');
speed_thresholded = speed<speed_threshold; % find events with speed below 5
if size(speed_thresholded,1)>size(speed_thresholded,2)
    speed_thresholded=speed_thresholded';
end

if size(speed_thresholded,1)<size(speed_thresholded,2)
    speed_thresholded=speed_thresholded';
end

if size(zscored_ripple,1)<size(zscored_ripple,2)
    zscored_ripple=zscored_ripple';
end

% speed_thresholded
thresholded = zscored_ripple > lowThresholdFactor;
thresholded = thresholded + speed_thresholded;
thresholded = thresholded-1;
thresholded(thresholded ~= 1) = 0;

% thresholded = zscored_ripple > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1
    start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1)
    stop(1) = [];
    start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass)
    ripples.onset =[];
    ripples.offset =[];
    ripples.peak_zscore =[];
    ripples.peaktimes = [];
    ripples.detectorinfo = [];
    disp('Detection by thresholding failed');
    return
else
    disp(['Detecting: ' num2str(length(firstPass)) ' ripple events after thresholding.']);
end

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*frequency;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
    if firstPass(i,1) - ripple(2) < minInterRippleSamples
        % Merge
        ripple = [ripple(1) firstPass(i,2)];
    else
        secondPass = [secondPass ; ripple];
        ripple = firstPass(i,:);
    end
end

secondPass = [secondPass ; ripple];
if isempty(secondPass)
    ripples.onset =[];
    ripples.offset =[];
    ripples.peak_zscore =[];
    ripples.peaktimes = [];
    ripples.detectorinfo = [];
    disp('Ripple merge failed');
    return
else
    disp(['Detecting: ' num2str(length(secondPass)) ' ripple events after merging.']);
end

% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
    [maxValue,maxIndex] = max(zscored_ripple([secondPass(i,1):secondPass(i,2)]));
    if maxValue > highThresholdFactor
        thirdPass = [thirdPass ; secondPass(i,:)];
        peakNormalizedPower = [peakNormalizedPower ; maxValue];
    end
end
if isempty(thirdPass),
    disp('Peak thresholding failed.');
    ripples.onset = [];
    ripples.offset = [];
    ripples.peaktimes = [];            %peaktimes? could also do these as timestamps and then ripples.ints for start/stops?
    ripples.peak_zscore = [];
    ripples.detectorinfo = [];
    return
else
    disp(['Detecting: ' num2str(length(thirdPass)) ' ripple events after ripple peak thresholding.']);
end

% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
    [minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
    peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

% Discard ripples that are way too long
ripples = [timevec(thirdPass(:,1)) timevec(peakPosition) ...
    timevec(thirdPass(:,2)) peakNormalizedPower];
duration = ripples(:,3)-ripples(:,1);
ripples(duration>maxRippleDuration/1000,:) = NaN;


% Discard ripples that are too short
ripples(duration<minRippleDuration/1000,:) = NaN;
ripples = ripples((all((~isnan(ripples)),2)),:);
disp(['Detecting: ' num2str(size(ripples,1)) ' ripple events.']);

% If a noise channel was provided, find ripple-like events and exclude them
bad = [];
if ~isempty(noise)

    noiselfp = filtfilt(b_ripple,1,noise);
    zscored_noise = zscore(abs(hilbert(noiselfp)));

    %     squaredNoise = bz_Filter(double(noise),'filter','butter','passband',passband,'order', 3).^2;
    % 	window = ones(windowLength,1)/windowLength;
    % 	zscored_noise = unity(Filter0(window,sum(squaredNoise,2)),sd,[]);

    excluded = logical(zeros(size(ripples,1),1));
    % Exclude ripples when concomittent noise crosses high detection threshold
    for i = 1:size(ripples,1)

        time_index = find(timevec > ripples(i,1) & timevec < ripples(i,3));

        if any(zscored_noise(time_index(1):time_index(2))>highThresholdFactor)
            excluded(i) = 1;
        end
    end
    bad = ripples(excluded,:);
    ripples = ripples(~excluded,:);
    disp(['Detecting: ' num2str(size(ripples,1)) ' ripple events after ripple-band noise removal.']);
end


% Optionally, plot results
if strcmp(show,'on')
    if plotType == 1
        figure('Name','Ripple-band filtered LFP');
        if ~isempty(noise)

            subplot(4,1,1)
            plot(timevec,signal)
            title('Ripple filtered signal LFP (125-300Hz)')
            subplot(4,1,2)
            plot(timevec,zscored_ripple)
            title('z-scored Ripple power (125-300Hz)')
            subplot(4,1,3)
            plot(timevec,noise)
            title('Ripple filtered noise LFP (125-300Hz)')
            subplot(4,1,4)
            plot(timevec,zscored_noise)
            title('z-scored Ripple power noise (125-300Hz)')
        else
            subplot(3,1,1)
            plot(timevec,signal)
            title('Ripple filtered signal LFP (125-300Hz)')
            xlabel('time (sec)')
            subplot(3,1,2)
            plot(timevec,zscored_ripple)
            title('z-scored Ripple power (125-300Hz)')
            xlabel('time (sec)')
            %             set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
            subplot(3,1,3)
            plot(timevec,speed)
            title('Movement (pixel change zscored)')
            xlabel('time (sec)')
            set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
        end

        figure('Name','Ripple detection');
        ripple_first = round(size(ripples,1)/2)-3;
        ripple_last = round(size(ripples,1)/2)+3;

        if ripple_last > size(ripples,1)
            ripple_last = size(ripples,1);
        end

        if ~isempty(ripples)
            if size(ripples,1)>10

                subplot(3,1,1)
                time_index = find(timevec>=ripples(ripple_first,1)-0.1 &timevec<=ripples(ripple_last,1)+0.1);
                plot(timevec(time_index),signal(time_index))
                title('Ripple filtered LFP (125-300Hz)')
                xlabel('time (sec)')
                set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

                subplot(3,1,2)
                plot(timevec(time_index),zscored_ripple(time_index))
                hold on
                plot([timevec(time_index(1)) timevec(time_index(end))],[highThresholdFactor highThresholdFactor],'k','linestyle','--');
                plot([timevec(time_index(1)) timevec(time_index(end))],[lowThresholdFactor lowThresholdFactor],'k');
                title('z-scored Ripple power (125-300Hz)')

                for j=ripple_first:ripple_last
                    hold on
                    plot([ripples(j,1) ripples(j,1)],[min(zscored_ripple(time_index)) ripples(j,4)],'g-');
                    plot([ripples(j,2) ripples(j,2)],[ripples(j,4) ripples(j,4)],'k-');
                    plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
                    plot([ripples(j,3) ripples(j,3)],[min(zscored_ripple(time_index)) ripples(j,4)],'r-');
                end
                hold on
                plot([timevec(time_index(1)) timevec(time_index(end))],[highThresholdFactor highThresholdFactor],'k','linestyle','--');
                plot([timevec(time_index(1)) timevec(time_index(end))],[lowThresholdFactor lowThresholdFactor],'k');
                xlabel('time (sec)')
                set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

                subplot(3,1,3)
                plot(timevec(time_index),speed(time_index))
                title('Movement (pixel change zscored)')
                hold on
                %         plot([timevec(time_index(1)) timevec(time_index(end))],[highThresholdFactor highThresholdFactor],'k','linestyle','--');
                for j=ripple_first:ripple_last
                    hold on
                    plot([ripples(j,1) ripples(j,1)],[min(speed(time_index)) lowThresholdFactor],'g-');
                    plot([ripples(j,2) ripples(j,2)],[lowThresholdFactor lowThresholdFactor],'k-');
                    plot([ripples(j,1) ripples(j,3)],[lowThresholdFactor lowThresholdFactor],'k-');
                    plot([ripples(j,3) ripples(j,3)],[min(speed(time_index)) lowThresholdFactor],'r-');
                end
                xlabel('time (sec)')
                set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
            end
        end
        % Check ripples
        %         fig = figure;
        %
        %         filter_type  = 'bandpass';
        %         passband = [0.5 3];
        %         filter_order = round(6*frequency/(max(passband)-min(passband)));  % creates filter for ripple
        %         norm_freq_range = passband/(frequency/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        %         b_slow = fir1(filter_order, norm_freq_range,filter_type);
        %         slow_signal = filtfilt(b_slow,1,lfp);
        %         zscored_ripple = zscore(abs(hilbert(slow_signal)));
        %         time_index = find(timevec>=ripples(ripple_first,1)-0.5 & timevec<=ripples(ripple_first,1)+0.5);
        %         plot(timevec(time_index),slow_signal(time_index),'k');hold on
        %         plot(timevec(time_index),signal(time_index),'r')
        %         title('LFP signal')
        %         title('Ripple filtered signal LFP (125-300Hz)')
        %         subplot(4,1,2)
        %         plot(timevec,zscored_ripple)
        %         title('z-scored Ripple power (125-300Hz)')
    end
end


%% BUZCODE Struct Output
temp = ripples; clear ripples

if ~isempty(temp)
    ripples.onset = temp(:,1);
    ripples.offset = temp(:,3);
    ripples.peaktimes = temp(:,2);            %peaktimes? could also do these as timestamps and then ripples.ints for start/stops?
    ripples.peak_zscore = temp(:,4);  %amplitudes?
else
    ripples.onset = [];
    ripples.offset = [];
    ripples.peaktimes = [];            %peaktimes? could also do these as timestamps and then ripples.ints for start/stops?
    ripples.peak_zscore = [];  %amplitudes?
end

ripples.best_channel = best_channel;  % save best channel info

% ripples.stdev = sd;
% if ~isempty(bad)
%     ripples.noise.times = bad(:,[1 3]);
%     ripples.noise.peaks = bad(:,[2]);
%     ripples.noise.peakNormedPower = bad(:,[4]);
% else
%     ripples.noise.times = [];
%     ripples.noise.peaks = [];
%     ripples.noise.peakNormedPower = [];
% end

%The detectorinto substructure
detectorinfo.detectorname = 'bz_FindRipples';
detectorinfo.detectiondate = datetime('today');
detectorinfo.detectionintervals = restrict;
detectorinfo.detectionparms = p.Results;
detectorinfo.detectionparms = rmfield(detectorinfo.detectionparms,'noise');

%Put it into the ripples structure
ripples.detectorinfo = detectorinfo;

%Save
if p.Results.saveMat
    save(fullfile(savepath,'extracted_ripples_events.mat'),'ripples')
end



%
% function y = Filter0(b,x)
%
% if size(x,1) == 1
% 	x = x(:);
% end
%
% if mod(length(b),2)~=1
% 	error('filter order should be odd');
% end
%
% shift = (length(b)-1)/2;
%
% [y0 z] = filter(b,1,x);
%
% y = [y0(shift+1:end,:) ; z(1:shift,:)];
%
% function [U,stdA] = unity(A,sd,restrict)
%
% if ~isempty(restrict),
% 	meanA = mean(A(restrict));
% 	stdA = std(A(restrict));
% else
% 	meanA = mean(A);
% 	stdA = std(A);
% end
% if ~isempty(sd),
% 	stdA = sd;
% end
%
% U = (A - meanA)/stdA;

