function units = load_session_clusters(paths, layout, cfg)
% LOAD_SESSION_CLUSTERS  Load Ellie's extracted_clusters_<sorter>.mat and label units.
%
% extracted_clusters_ks4.mat -> clusters_ks4 (one probe):
%   .cluster_id    1-based (SpikeInterface id + 1)
%   .spike_times   s, SpikeGLX time of this recording; .spike_id = cluster_id per spike
%   .peak_channel  channel-map row of the peak waveform
%   quality metrics (amplitude_cutoff, sliding_rp_violation, snr, ...)
%
% units adds, per cluster:
%   .shank / .depth  from the channel map at peak_channel
%   .region          "V1" (L5 + cfg.depth_ranges.V1), "CA1" (CA1 + cfg.depth_ranges.CA1)
%                    or "other", using that shank's PSD depths
%   .is_good         passes create_cluster_selection_params('sorting_option','ellie')

f = fullfile(paths.analysis_dir, sprintf('extracted_clusters_%s.mat', cfg.sorter));
S = load(f);
C = S.(['clusters_' cfg.sorter]);
C = C(1);

[~, good] = select_clusters(C, create_cluster_selection_params('sorting_option', 'ellie'));

units.cluster_id   = C.cluster_id(:);
units.shank        = layout.kcoords(C.peak_channel(:));
units.depth        = layout.ycoords(C.peak_channel(:));
units.peak_channel = C.peak_channel(:);
units.firing_rate  = C.firing_rate(:);
units.cell_type    = C.cell_type(:);
units.is_good      = logical(good(:));
units.region       = repmat("other", numel(units.cluster_id), 1);
for s = layout.shank
    on_shank = units.shank == s.id;
    units.region(on_shank & units.depth >= s.CA1_range(1) & units.depth <= s.CA1_range(2)) = "CA1";
    units.region(on_shank & units.depth >= s.V1_range(1)  & units.depth <= s.V1_range(2))  = "V1";
end
units.spike_times = C.spike_times(:);
units.spike_id    = C.spike_id(:);
units.source_file = f;

fprintf('Units: %d total | V1 %d (%d good) | CA1 %d (%d good)\n', numel(units.cluster_id), ...
    sum(units.region == "V1"), sum(units.region == "V1" & units.is_good), ...
    sum(units.region == "CA1"), sum(units.region == "CA1" & units.is_good));
end
