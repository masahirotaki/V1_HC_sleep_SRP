%% Stimulus decoding: what CA1 and V1 population spiking carries about the GAVNIK ABCD sequence
% Reads Ellie's data from Y:\Ellie\DATA and writes into each recording's own folder,
% Y:\Ellie\DATA\SUBJECTS\<subject>\analysis\<date>\<stimulus>\Decoding\ (decoding_* files).
%
% Steps (each can be switched off below):
%   1. decoding        run_abcd_decoding: element identity (A/B/C/D and A/B/C/D/grey) and time
%                      within the sequence; decoders: Poisson Bayes, multinomial Bayes
%                      (count-conditioned, gain-invariant), shrinkage LDA, PCA+LDA, PLS-DA,
%                      linear SVM; 10-fold CV over sequences, label-shuffle nulls; single-unit
%                      time-cell plots, PCA/PLS. Also saves decoder templates for sleep.   ~16 min
%   2. multinomial     plot_multinomial_results: posterior figures of the multinomial decoder
%   3. pca_fa          run_abcd_pca_fa: PCA and factor analysis (DataHigh FA, q by CV) single-
%                      trial and mean trajectories                                         ~10 min
%   4. seqlik          run_abcd_sequence_likelihood: Lseq (Huh et al. 2026) order matrix from
%                      training sequences, held-out sequences vs grey windows, 10-fold CV  ~8 min
%   5. seqlik_figures  example rasters (CA1) and the trial-structure plot
%   6. contributions   run_abcd_seqlik_contributions + per-spike view: which neurons carry Lseq ~5 min
%   7. assembly_cv     run_abcd_assembly_cv: ICA assemblies (20 ms bins) + LDA, 10-fold CV of
%                      A/B/C/D/grey, the awake test of the sleep reactivation readout
% Steps 5-6 need step 4; steps 2 needs step 1.
%
% Defaults: 50 ms time bins, 200 permutations, 1000 Lseq shuffles, original 'random' grey
% windows; A/B/C/D identity counted from each element's onset to its offset ([0, 150 ms] for
% GAVNIK; until 2026-10-15 it was [30, 180] ms, so earlier saved results used that window).
% Windows follow each session's measured timing (get_abcd_trials: element on-time, onset
% spacing, sequence length), so GAVNIK250 and TRAIN sessions run with the same code.
% Spike sorting is one sort per day, so cluster ids match the Sleep_Box recordings.

clear; close all
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
cfg = project_config();

%% Recordings
% One recording:
sessions = struct('subject', {'M00069'}, 'date', {'20250929'}, 'stimulus', {'GAVNIK_ABCD_1'});
% Every animal's first novel GAVNIK / GAVNIK250 day (config/first_novel_sessions.m). The
% TRAIN day (M00013: 150 ms gratings every 300 ms) also runs (get_abcd_trials reads each
% session's timing) but is left out for consistency; drop the contains() line to include it.
S = first_novel_sessions();
S = S(contains({S.novel_stimulus}, 'GAVNIK'));
sessions = struct('subject', {S.subject}, 'date', {S.date}, 'stimulus', {S.novel_stimulus});

%% Steps
steps.decoding      = true;
steps.multinomial   = true;
steps.pca_fa        = true;
steps.seqlik        = true;
steps.seqlik_figures = true;
steps.contributions = true;
steps.assembly_cv   = true;

%% Options (empty struct = defaults)
opts_decoding = struct();   % e.g. struct('time_bin', 0.075) or struct('n_perm', 500)
opts_pca_fa   = struct();   % e.g. struct('bin', 0.02, 'sigma', 0.04)
opts_seqlik   = struct();   % e.g. struct('grey_sampling', 'nonoverlap') for non-overlapping grey windows
opts_contrib  = struct();   % e.g. struct('n_shuffle', 1000)
opts_assembly = struct();   % e.g. struct('bin', 0.02)

%% Run
status = table('Size', [0 6], 'VariableTypes', {'string', 'string', 'string', 'string', 'logical', 'string'}, ...
    'VariableNames', {'subject', 'date', 'stimulus', 'step', 'ok', 'message'});

for i = 1:numel(sessions)
    s = sessions(i);
    fprintf('\n################ %s %s %s ################\n', s.subject, s.date, s.stimulus);
    decoding_file = fullfile(cfg.subjects_root, s.subject, 'analysis', s.date, s.stimulus, 'Decoding', ...
        sprintf('decoding_%s_results.mat', s.stimulus));
    todo = {
        'decoding',       steps.decoding,       @() run_abcd_decoding(s.subject, s.date, s.stimulus, cfg, opts_decoding)
        'multinomial',    steps.multinomial,    @() plot_multinomial_results(decoding_file, cfg)
        'pca_fa',         steps.pca_fa,         @() run_abcd_pca_fa(s.subject, s.date, s.stimulus, cfg, opts_pca_fa)
        'seqlik',         steps.seqlik,         @() run_abcd_sequence_likelihood(s.subject, s.date, s.stimulus, cfg, opts_seqlik)
        'seqlik_figures', steps.seqlik_figures, @() seqlik_figures(s, cfg)
        'contributions',  steps.contributions,  @() contributions(s, cfg, opts_contrib)
        'assembly_cv',    steps.assembly_cv,    @() run_abcd_assembly_cv(s.subject, s.date, s.stimulus, cfg, opts_assembly)
        };
    for k = 1:size(todo, 1)
        if ~todo{k, 2}; continue; end
        status(end + 1, :) = run_step(s, todo{k, 1}, todo{k, 3}); %#ok<SAGROW>
        close all
    end
end
disp(status);


%% ------------------------------------------------------------------------
function row = run_step(s, name, fn)
t0 = tic;
fprintf('\n---- %s ----\n', name);
try
    fn();
    row = {string(s.subject), string(s.date), string(s.stimulus), string(name), true, ...
        string(sprintf('%.1f min', toc(t0) / 60))};
catch err
    fprintf(2, '%s %s %s: %s failed: %s\n', s.subject, s.date, s.stimulus, name, err.message);
    row = {string(s.subject), string(s.date), string(s.stimulus), string(name), false, string(err.message)};
end
end

function seqlik_figures(s, cfg)
plot_seqlik_examples(s.subject, s.date, s.stimulus, 'CA1_all', cfg);
plot_seqlik_examples(s.subject, s.date, s.stimulus, 'CA1_good', cfg);
plot_seqlik_trial_structure(s.subject, s.date, s.stimulus, cfg);
end

function contributions(s, cfg, opts)
run_abcd_seqlik_contributions(s.subject, s.date, s.stimulus, cfg, opts);
plot_seqlik_contrib_per_spike(s.subject, s.date, s.stimulus, cfg);
end
