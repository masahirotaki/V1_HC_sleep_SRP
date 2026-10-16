function [selected_clusters,cluster_filter_index] = select_clusters(clusters,params)
if nargin < 2
    % default values of the params
    params.isi_violations_ratio = @(x) x<=0.1;
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude_median = @(x) x>50; %IBL 50
    params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
    params.drift_ptp = @(x) x<5;
    params.drift_std = @(x) x<1;
    params.drift_mad = @(x) x<1;
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
    params.amplitude_cv_median = @(x) x<=0.7;
    params.amplitude_cv_range = @(x) x<=0.7;
end
selected_clusters = struct();

excluded_clusters_by_metric = struct(); % 1/3 initialise log of filtered out clusters
all_metrics = fieldnames(params);
cluster_filter_index = ones(size(clusters.cluster_id));
for iMetric = 1:size(all_metrics,1)
    temp_cluster_index = params.(all_metrics{iMetric})(clusters.(all_metrics{iMetric}));
    if ~isstring(clusters.(all_metrics{iMetric})) % if ischar (region label) skip this
        if sum(isnan(clusters.(all_metrics{iMetric}))) > 0
            temp_cluster_index(isnan(clusters.(all_metrics{iMetric}))) = 1;
        end
    end
    disp([num2str(sum(temp_cluster_index)), ' clusters satisfy threshold of ',all_metrics{iMetric}])

    failed = ~temp_cluster_index; % 2a/3 store cluster IDs that failed this metric
    excluded_clusters_by_metric.(all_metrics{iMetric}) = clusters.cluster_id(failed); % 2b/3 store cluster IDs that failed this metric
    cluster_filter_index = cluster_filter_index &temp_cluster_index;
end
disp(['Overall ',num2str(sum(cluster_filter_index)), ' clusters satisfy all metrics specified'])
all_cluster_fields = fieldnames(clusters);

for iField = 1:size(all_cluster_fields,1)
    if contains(all_cluster_fields{iField},'timevec')
        selected_clusters.timevec = clusters.timevec;
    elseif contains(all_cluster_fields{iField},'merged_spike_id')
        spike_id_selected = ismember(clusters.spike_id,clusters.cluster_id(cluster_filter_index));
        selected_clusters.merged_spike_id = clusters.merged_spike_id(spike_id_selected);
    elseif strcmp(all_cluster_fields{iField},'spike_id')
        spike_id_selected = ismember(clusters.spike_id,clusters.cluster_id(cluster_filter_index));
        selected_clusters.spike_id = clusters.spike_id(spike_id_selected);
        %     elseif contains(all_cluster_fields{iField},'unstable_ids')
        %         selected_clusters.unstable_ids = clusters.unstable_ids;
    elseif contains(all_cluster_fields{iField},'spike_times')
        selected_clusters.spike_times = clusters.spike_times(spike_id_selected);
    elseif contains(all_cluster_fields{iField},'probe_id')
        selected_clusters.probe_id = clusters.probe_id;
    elseif contains(all_cluster_fields{iField},'sorter')
        selected_clusters.sorter = clusters.sorter;
    elseif contains(all_cluster_fields{iField},'probe_hemisphere')
        selected_clusters.probe_hemisphere = clusters.probe_hemisphere;
    elseif contains(all_cluster_fields{iField},'lap_seesion_ID')
        selected_clusters.lap_seesion_ID = clusters.lap_seesion_ID;
    elseif contains(all_cluster_fields{iField},'pass_landmarks_laps')
        selected_clusters.pass_landmarks_laps = clusters.pass_landmarks_laps;
    elseif contains(all_cluster_fields{iField},'spatial_response_all_2cm')
        selected_clusters.spatial_response_all_2cm = clusters.spatial_response_all_2cm(cluster_filter_index,:,:);
    elseif contains(all_cluster_fields{iField},'excluded_by_metric')
        % selected_clusters.excluded_by_metric = clusters.excluded_by_metric;

    elseif contains(all_cluster_fields{iField},'params')

    elseif iscell(clusters.(all_cluster_fields{iField})) % if it is a cell structure (Usually for RF which is cell structure)

        if size(clusters.(all_cluster_fields{iField}),2)>1 & size(clusters.(all_cluster_fields{iField}),1)>1
            temp_cluster_field = clusters.(all_cluster_fields{iField});
            selected_clusters.(all_cluster_fields{iField}) = {temp_cluster_field{cluster_filter_index,:}};
            selected_clusters.(all_cluster_fields{iField})=reshape(selected_clusters.(all_cluster_fields{iField}),[sum(cluster_filter_index) size(temp_cluster_field,2)]);
        else
            temp_cluster_field = clusters.(all_cluster_fields{iField});
            if length(temp_cluster_field)==length(cluster_filter_index)
                selected_clusters.(all_cluster_fields{iField}) = {temp_cluster_field{cluster_filter_index}}';
            else
                selected_clusters.(all_cluster_fields{iField})=temp_cluster_field;
            end
        end
    else
        temp_cluster_field = clusters.(all_cluster_fields{iField});
        selected_clusters.(all_cluster_fields{iField}) = temp_cluster_field(cluster_filter_index,:);
    end
end

selected_clusters.excluded_by_metric = excluded_clusters_by_metric; % 3/3
selected_clusters.params = params;% Save the parameters that were used for selection.