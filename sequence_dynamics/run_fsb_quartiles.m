% RUN_FSB_QUARTILES  Broad V1 and CA1 units by quartile of rate-corrected spontaneous FSB.
%
% Re-plots from the saved per-day results (run_burst_isi5_fr5 / run_burst_isi5_fr5_CA1
% must have been run) into <root>\fsb_quartiles\<smoothing>\ (see plot_fsb_quartiles), with
% the per-unit quartile table as fsb_quartile_units.csv and the per-unit FSB distributions
% (plot_fsb_distribution) in <root>\fsb_quartiles\. Reads only those data folders.

roots = {'Y:\Ellie\DATA\visual_sequence_oscillation\burst_isi5ms_fr5Hz', 'burst_psth_M*.mat'
         'Y:\Ellie\DATA\visual_sequence_oscillation\burst_isi5ms_fr5Hz_CA1', 'burst_psth_CA1_*.mat'};
smoothing = {'smoothed_10ms', 0.010; 'unsmoothed', 0};
min_exp = 5;   % minimum expected Poisson burst spikes per unit

setup_paths();
for k = 1:size(roots, 1)
    F = dir(fullfile(roots{k, 1}, 'data', roots{k, 2}));
    R = cell(numel(F), 1);
    for i = 1:numel(F)
        tmp = load(fullfile(F(i).folder, F(i).name));
        R{i} = tmp.r;
    end
    R = [R{:}];
    out = fullfile(roots{k, 1}, 'fsb_quartiles');
    for s = 1:size(smoothing, 1)
        Q = plot_fsb_quartiles(R, fullfile(out, smoothing{s, 1}), smoothing{s, 2}, min_exp);
    end
    writetable(Q, fullfile(out, 'fsb_quartile_units.csv'));
    plot_fsb_distribution(R, Q, out, min_exp);   % figQ_fsb_distribution (independent of smoothing)
end
