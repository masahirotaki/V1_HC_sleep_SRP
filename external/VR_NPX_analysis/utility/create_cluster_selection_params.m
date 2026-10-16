function params = create_cluster_selection_params(varargin)

% Default values
p = inputParser;
addParameter(p,'sorting_option','spikeinterface',@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))

% assign parameters (either defaults or given)
parse(p,varargin{:});
sorting_option = p.Results.sorting_option;

if contains(sorting_option,'spikeinterface')
    % default values of the params
%     params.isi_violations_ratio = @(x) x<=0.1;
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude_median = @(x) abs(x)>50; %IBL 50
    params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
    % params.drift_ptp = @(x) x<5;
    % params.drift_std = @(x) x<1;
    % params.drift_mad = @(x) x<1;
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
    params.amplitude_cv_median = @(x) x<=0.7;
    params.amplitude_cv_range = @(x) x<=0.7;

elseif contains(sorting_option,'old')
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude = @(x) x>50; %IBL 50
    params.isi_viol = @(x) x<=0.1;
    params.duration = @(x) x<=0.8 & x>= 0.2;
elseif contains(sorting_option,'masa') 
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude_median = @(x) abs(x)>20; %IBL 50 but bombcell 20
    params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0001; % Peters et al 2022
    params.amplitude_cv_median = @(x) x<=0.7;
    params.amplitude_cv_range = @(x) x<=0.7;
    params.firing_rate = @(x) x>= 0.01;% Due to sleep/wave difference, we
    %can have silent neuron during sleep or awake
    % params.presence_ratio = @(x) x>= 0.3;% Due to sleep/wave difference,
    %we can have silent neuron during sleep or awake during sleep or awake
    params.sd_ratio = @(x) x<=4;
    params.snr = @(x) x>= 2;
    params.half_width = @(x) x<= 0.0003;
elseif contains(sorting_option,'ellie')
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units. A histogram of a unit's spike amplitudes should have a Gaussian distribution. 
    % Deviations from this distribution are used to estimate the number of spikes missing from the unit (false negative rate). The upper tail of the distribution is used to estimate the fraction of spikes missing from the lower tail). Can be computed on chunks for a long recording as drift can impact amplitudes and hence the normalness of the distribution. 
    params.amplitude_median = @(x) abs(x)>20; %IBL 50 but bombcell 20. Larger value = larger signal, better unit.
    params.sliding_rp_violation = @(x) x<=0.1; % IBL uses 10%. Refractory period violations, similar to ISI violations but the sliding window allows inclusion of cells with burst firing (useful for HPC); 0.1 = 10% of the expected number of violations under chance firing.
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0001; % Peters et al 2022 % the duration between negative and positive peaks
    params.amplitude_cv_median = @(x) x<=0.7; % cv = coefficient of variation, variability around the mean; ratio of s.d./mean amplitude. A low median value indicates a stereotypical spike shape i.e. a well-isolated unit.
    params.amplitude_cv_range = @(x) x<=0.7; % cv range may be high if there is noise contamination
    params.firing_rate = @(x) x>= 0.01;% number of spikes divided by number of seconds in the recording. Due to sleep/wake difference, we
    %can include neurons that are silent during sleep or wake
    params.presence_ratio = @(x) x>= 0.3;% Due to sleep/wake difference, we can include neurons that are silent during sleep or wake. 
    % if you are analysing a short segment of the recording, you could disable the default presence ratio filter, to maximise the number of units available.
    % if it is not clear whether a unit has a low presence ratio due to drift or selective firing, plot spike amplitude over time; if amplitude drifts around the spike detection threshold, there can be a very variable firing rate
    params.sd_ratio = @(x) x<=4; % at the peak of the spike, the sd of the voltage should be that of noise (in which case sd_ratio = 1); contamination will increase the sd.
    params.snr = @(x) x>= 2; % Ratio of the max amplitude of the mean spike waveform to the sd of the background noise on one (the peak) channel.
    params.half_width = @(x) x<= 0.0003; %halfwidth is the duration in seconds at 50% amplitude

else

end

end