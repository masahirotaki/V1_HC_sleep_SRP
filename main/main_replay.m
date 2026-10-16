%% Replay of the ABCD response order (sequence likelihood, Huh et al. 2026) in CA1 and V1
% Needs main_stimulus_decoding (step 'seqlik') for the same recordings: the order matrices,
% folds and units come from decoding_<stimulus>_seqlik_results.mat.
%
% Block 1  Awake replay during the stimulus recording (run_abcd_awake_replay)
%          - awake ripples (FindRipples_masa, cfg.ripple) and MUA candidate events
%            (detect_candidate_events_masa: CA1 MUA z >= 3, edges at 0, >= 100 ms, merged
%            < 50 ms, >= 5 neurons) with ripple-band peak z >= 5
%          - CA1 spikes scored with the CA1 order, V1 spikes (same events) with the V1 order;
%            cell-ID shuffle (alpha 0.025 per tail) + matched-duration grey calibration;
%            events overlapping a training sequence use the CV fold that left it out
%          - timing of ripples / candidates / replay around A onset (-0.5 to 1 s)
%          - number of ripples inside each candidate; ripples inside significant events
%          - peri-event population PSTHs (CA1 and V1) for all candidates, CA1-significant and
%            V1-significant events, all and grey-screen-only, vs matched grey times
%          - stimulus-element assembly reactivation in grey-screen awake ripples
%            (run_assembly_reactivation, see Block 2)
%          Outputs: <stimulus>\replay\replay_<stimulus>_*            (~2 min per recording)
%
% Block 2  Sleep (set run_sleep = true): ripple clusters + MUA bursts, Lseq replay of ripples /
%          clusters / MUA bursts (CA1 and V1), ripple timing around UP / DOWN states,
%          peri-event LFP / CSD / spiking, and assembly reactivation. Needs
%          main_sleep_preprocessing for each Sleep_Box.
%          Outputs: <Sleep_Box_x>\replay\replay_<Sleep_Box_x>_*  and
%                   <Sleep_Box_x>\detection\extracted_peri_event_signals.mat / _csd.pdf

clear; close all
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
cfg = project_config();

%% Recordings
% One day: the awake stimulus recording and the sleep recordings before / after it
% sessions = struct('subject', {'M00069'}, 'date', {'20250929'}, 'stimulus', {'GAVNIK_ABCD_1'}, ...
%     'pre_sleep', {{'Sleep_Box_1'}}, 'post_sleep', {{'Sleep_Box_2'}});
% Every animal's first novel GAVNIK day:
S = first_novel_sessions();
S = S(contains({S.novel_stimulus}, 'GAVNIK'));   % GAVNIK days only (drop this line to add TRAIN)
sessions = struct('subject', {S.subject}, 'date', {S.date}, 'stimulus', {S.novel_stimulus}, ...
    'pre_sleep', {S.pre_sleep}, 'post_sleep', {S.post_sleep});

%% Block 1: awake replay (CA1 and V1)
opts_awake = struct();   % e.g. struct('alpha', 0.05) for the paper's 10% two-tailed convention,
                         %      struct('psth_significance', 'calibrated'), struct('psth_window', [-1 1])
opts_assembly = struct();   % run_assembly_reactivation options, e.g. struct('n_shuffle', 2000)
run_awake = true;
if run_awake
    for i = 1:numel(sessions)
        s = sessions(i);
        fprintf('\n######## awake replay: %s %s %s ########\n', s.subject, s.date, s.stimulus);
        try
            run_abcd_awake_replay(s.subject, s.date, s.stimulus, cfg, opts_awake);
        catch err
            fprintf(2, '%s %s %s awake replay failed: %s\n', s.subject, s.date, s.stimulus, err.message);
        end
        close all
        try
            run_assembly_reactivation(s.subject, s.date, s.stimulus, s.stimulus, cfg, opts_assembly);
        catch err
            fprintf(2, '%s %s %s awake assembly reactivation failed: %s\n', s.subject, s.date, s.stimulus, err.message);
        end
        close all
        try
            analyze_replay_cosignificance(s.subject, s.date, s.stimulus, s.stimulus, cfg);
        catch err
            fprintf(2, '%s %s %s awake co-significance failed: %s\n', s.subject, s.date, s.stimulus, err.message);
        end
        close all
    end
end

%% Block 2: sleep (pre and post): clusters, replay, UP/DOWN timing, peri-event signals
% Needs main_sleep_preprocessing (run_sleep_detection) on every Sleep_Box first.
% Per sleep recording:
%   analyze_ripple_mua_overlap  ripple clusters (Voroslakos et al. 2026) and CA1 MUA bursts
%                               (1 ms, sigma 5 ms, peak z > 3, edges z = 1, <= 300 ms)
%   run_sleep_replay            Lseq (cell-ID shuffle, alpha 0.025 per tail + matched controls)
%                               for each ripple [onset, max(offset, onset + 100 ms)], each merged
%                               ripple cluster and MUA bursts without ripple; CA1 and V1 templates
%                               from the awake stimulus recording; continuation tests in clusters
%   analyze_ripple_updown       ripples around UP / DOWN onsets and at normalised positions in UP
%                               (<= 2 s) and DOWN states; split by solo / cluster / replay
%   compute_peri_event_signals  peri-event LFP / CSD per shank and CA1 / V1 spiking around ripples,
%                               UP onsets and DOWN onsets (saved to detection\ for later analyses)
%   run_assembly_reactivation   awake stimulus-element ICA assemblies (20 ms, z per neuron + common-
%                               mode removal); pass 1: per ripple, max over bins of total activation
%                               vs whole-ripple cell-ID shuffle; pass 2: element (A-D/grey) per bin
%                               in significant ripples (LDA on assemblies, multinomial on counts)
%   analyze_replay_cosignificance  % of events with significant Lseq in CA1 / V1, CA1-V1 co-
%                               significance vs independence and vs a matched-size null, example
%                               rasters (also run on the awake candidates in Block 1)
opts_sleep = struct();   % run_sleep_replay options, e.g. struct('states', {{'SWS'}})
run_sleep = true;
if run_sleep
    for i = 1:numel(sessions)
        s = sessions(i);
        for sleep_name = [s.pre_sleep s.post_sleep]
            rec = sleep_name{1};
            fprintf('\n######## sleep: %s %s %s (template %s) ########\n', s.subject, s.date, rec, s.stimulus);
            steps = {
                @() analyze_ripple_mua_overlap(s.subject, s.date, rec, cfg)
                @() run_sleep_replay(s.subject, s.date, rec, s.stimulus, cfg, opts_sleep)
                @() analyze_ripple_updown(s.subject, s.date, rec, cfg, struct('stimulus', s.stimulus))
                @() compute_peri_event_signals(s.subject, s.date, rec, cfg)
                @() run_assembly_reactivation(s.subject, s.date, s.stimulus, rec, cfg, opts_assembly)
                @() analyze_replay_cosignificance(s.subject, s.date, rec, s.stimulus, cfg)
                };
            for k = 1:numel(steps)
                try
                    steps{k}();
                catch err
                    fprintf(2, '%s %s %s step %d failed: %s\n', s.subject, s.date, rec, k, err.message);
                end
                close all
            end
        end
    end

    % pre vs post: fraction of calibrated significant events per state and population
    for i = 1:numel(sessions)
        s = sessions(i);
        for sleep_name = [s.pre_sleep s.post_sleep]
            f = fullfile(cfg.subjects_root, s.subject, 'analysis', s.date, sleep_name{1}, 'replay', ...
                sprintf('replay_%s_template_%s_summary.csv', sleep_name{1}, s.stimulus));
            if ~isfile(f); continue; end
            T = readtable(f);
            fprintf('\n%s %s %s:\n', s.subject, s.date, sleep_name{1});
            disp(T);
        end
    end
end
