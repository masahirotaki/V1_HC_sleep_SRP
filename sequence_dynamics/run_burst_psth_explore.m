function run_burst_psth_explore(out_dir, recompute)
% RUN_BURST_PSTH_EXPLORE  Exploratory look at burst vs single-spike dynamics in V1
% during the day-5 test sequences (ABCD, A_CD, E_CD) of the 150 ms and 250 ms protocols.
%
% Reads Ellie's data only. Per-session results are cached in <out_dir>/cache and
% figures are written to <out_dir> (a local folder, not the server).

if nargin < 2; recompute = false; end
setup_paths();
cfg = project_config();
cache_dir = fullfile(out_dir, 'cache');
if ~exist(cache_dir, 'dir'); mkdir(cache_dir); end

S = day5_test_sessions();
R = cell(numel(S), 1);
for i = 1:numel(S)
    f = fullfile(cache_dir, sprintf('burst_psth_%s_%s_%d.mat', S(i).subject, S(i).date, S(i).protocol));
    if exist(f, 'file') && ~recompute
        tmp = load(f); R{i} = tmp.r;
        continue
    end
    try
        r = compute_session_burst_psth(S(i).subject, S(i).date, S(i).protocol, cfg);
        r.week = S(i).week;
        save(f, 'r', '-v7.3');
        R{i} = r;
    catch err
        fprintf(2, 'FAILED %s %s %d: %s\n', S(i).subject, S(i).date, S(i).protocol, err.message);
    end
end
R = [R{~cellfun(@isempty, R)}];

for p = [150 250]
    Rp = R([R.protocol] == p);
    if isempty(Rp); continue; end
    plot_burst_psth_explore(Rp, out_dir, sprintf('%dms_pooled', p));
    for i = 1:numel(Rp)
        plot_burst_psth_explore(Rp(i), fullfile(out_dir, 'per_session'), ...
            sprintf('%dms_%s_%s', p, Rp(i).subject, Rp(i).date), true);
    end
end
end
