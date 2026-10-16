function paths = get_session_paths(subject, date, stimulus_name, cfg)
% GET_SESSION_PATHS  Locate all input files for one recording and define outputs.
%
%   paths = get_session_paths('M00069', '20250929', 'Sleep_Box_1', cfg)
%
% Inputs (read-only, Ellie's server):
%   analysis_dir  SUBJECTS/<subj>/analysis/<date>/<stimulus>   extracted_behaviour/clusters
%   psd_file      SUBJECTS/<subj>/analysis/<date>/depths_from_PSD.mat      CA1 / L5 depth per shank
%   csd_file      SUBJECTS/<subj>/analysis/<date>/earliest_V1sink_CSD.mat  L4 depth per shank
%   imec_dir      SUBJECTS/<subj>/ephys/<date>/<date>_<sid>/<run>_g<gs>/<run>_g<gs>_imec<probe>
%   lf_bin_file   *tcat*.lf.bin (CatGT LFP, 2.5 kHz) or *.lf.bin
% Outputs (written into the recording's analysis folder, see check_output_path):
%   output_dir    = analysis_dir
%   detection_dir <output_dir>/detection

date = char(string(date));
stimulus_name = char(stimulus_name);
rec = find_recording_in_table(cfg.stim_table, subject, date, stimulus_name);

paths.subject   = subject;
paths.date      = date;
paths.stimulus  = stimulus_name;
paths.session_id = rec.session_id;
paths.gs_number  = rec.gs_number;

%% Ellie's processed data
date_analysis_dir = fullfile(cfg.subjects_root, subject, 'analysis', date);
paths.analysis_dir = fullfile(date_analysis_dir, stimulus_name);
assert(exist(paths.analysis_dir, 'dir') == 7, 'Analysis folder not found: %s', paths.analysis_dir);
paths.psd_file = fullfile(date_analysis_dir, 'depths_from_PSD.mat');
paths.csd_file = fullfile(date_analysis_dir, 'earliest_V1sink_CSD.mat');

%% Raw ephys
ephys_date_dir = fullfile(cfg.subjects_root, subject, 'ephys', date);
if strlength(rec.session_id) == 0
    sess_dir = ephys_date_dir;
else
    candidates = {fullfile(ephys_date_dir, [date '_' char(rec.session_id)]), ...
                  fullfile(ephys_date_dir, [subject '_' date '_' char(rec.session_id)])};
    found = find(cellfun(@(d) exist(d, 'dir') == 7, candidates), 1);
    assert(~isempty(found), 'Ephys session folder not found: %s', candidates{1});
    sess_dir = candidates{found};
end
d = dir(sess_dir);
is_g = [d.isdir] & ~cellfun(@isempty, regexp({d.name}, sprintf('_g%d$', rec.gs_number), 'once'));
assert(nnz(is_g) == 1, 'Expected one *_g%d folder in %s, found %d', rec.gs_number, sess_dir, nnz(is_g));
g_name = d(is_g).name;
paths.imec_dir = fullfile(sess_dir, g_name, sprintf('%s_imec%d', g_name, cfg.probe_id));

lf = dir(fullfile(paths.imec_dir, '*tcat*.lf.bin'));
if isempty(lf); lf = dir(fullfile(paths.imec_dir, '*.lf.bin')); end
assert(~isempty(lf), 'No lf.bin found in %s', paths.imec_dir);
paths.lf_bin_file  = fullfile(lf(1).folder, lf(1).name);
paths.lf_meta_file = regexprep(paths.lf_bin_file, '\.bin$', '.meta');

cm = dir(fullfile(paths.imec_dir, '*ChanMap*.mat'));
assert(~isempty(cm), 'No channel map found in %s', paths.imec_dir);
paths.chanmap_file = fullfile(cm(1).folder, cm(1).name);

%% Outputs
paths.output_dir    = paths.analysis_dir;
paths.detection_dir = fullfile(paths.output_dir, 'detection');
end


function rec = find_recording_in_table(table_file, subject, date, stimulus_name)
% Read the subject's sheet with readcell so mixed numeric/text session_id survives.
% A stimulus listed k times on one day is addressed as <name>_<k>, matching
% Ellie's analysis folder names (e.g. OP_Tuning_2).
C = readcell(table_file, 'Sheet', subject);
hdr = string(C(1, :));
body = C(2:end, :);
column = @(name) body(:, find(strcmpi(hdr, name), 1));

dates = cellfun(@cell_to_double, column('date'));
stims = cellfun(@cell_to_string, column('stimulus_type'));
same_day = dates == str2double(date);
row = find(same_day & stims == string(stimulus_name));
if numel(row) > 1
    error('%s is listed %d times on %s %s; use %s_<k> (k = 1..%d)', ...
        stimulus_name, numel(row), subject, date, stimulus_name, numel(row));
elseif isempty(row)
    tok = regexp(stimulus_name, '^(.*)_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        rows = find(same_day & stims == string(tok{1}));
        k = str2double(tok{2});
        if k <= numel(rows); row = rows(k); end
    end
end
assert(~isempty(row), '%s %s %s not found in %s', subject, date, stimulus_name, table_file);

sid = column('session_id');
gs  = column('gs_number');
rec.session_id = cell_to_string(sid{row});
rec.gs_number  = cell_to_double(gs{row});
end


function x = cell_to_double(c)
if isnumeric(c); x = double(c);
elseif ismissing(c); x = NaN;
else; x = str2double(c);
end
end

function s = cell_to_string(c)
if isnumeric(c) && ~isnan(c); s = string(num2str(c));
elseif (isnumeric(c) && isnan(c)) || any(ismissing(c)); s = "";
else; s = strtrim(string(c));
end
end
