function trials = get_abcd_trials(paths, element_names)
% GET_ABCD_TRIALS  Sequence and element timing of a 4-grating sequence recording:
% GAVNIK(250)_ABCD (gratings back to back) or TRAIN (150 ms gratings separated by 150 ms grey).
%
% extracted_task_info.mat -> Task_info:
%   .stim_onset        one entry per ELEMENT onset (s, SpikeGLX time)
%   .stim_orientation  grating orientation of that element (deg; radians in older sessions,
%                      e.g. M00013, converted here)
%   .pd_on/.pd_off.sglxTime  photodiode transitions. GAVNIK: the patch flips at every element
%                      (pd_on = A and C onsets, pd_off = B and D onsets). TRAIN: one on/off
%                      pulse per element (pd_off = element offset).
% Elements are grouped into sequences at onset gaps > opts gap (0.7 s); sequences without
% exactly 4 elements (e.g. a missed onset) are dropped with a warning.
% Element offsets: the photodiode pd_off after each onset when there is one pulse per element;
% otherwise (back to back) the next element's onset, and for the last element its onset + the
% median onset-to-onset interval.
%
% trials (one row per sequence):
%   .element_onset / .element_offset  [nSeq x nElem] s
%   .seq_onset / .seq_offset          [nSeq x 1] s (A onset / D offset)
%   .orientation                      [1 x nElem] deg, identical for every sequence
%   .element_duration                 median on-time of an element (s)
%   .soa_median                       median onset-to-onset interval within a sequence (s)
%   .element_start / .element_end     [1 x nElem] nominal element on/off re A onset
%                                     ((e-1) * soa_median, + element_duration)
%   .seq_duration                     nominal A onset -> D offset (s)
%   .isi_median                       grey between elements (0 for back-to-back designs)
%   .block                            presentation block (pauses > 5 s split blocks)

if nargin < 2; element_names = {'A', 'B', 'C', 'D'}; end
n_elem = numel(element_names);
gap_new_seq = 0.7;   % s between consecutive onsets that starts a new sequence
T = load(fullfile(paths.analysis_dir, 'extracted_task_info.mat'));
T = T.Task_info;
on = T.stim_onset(:);
ori = T.stim_orientation(:);
if all(ori >= 0 & ori <= 2 * pi + 1e-6) && any(abs(ori - round(ori)) > 1e-6)
    ori = round(rad2deg(ori), 4);          % radians in older sessions
end

% group into sequences
seq_id = cumsum([1; diff(on) > gap_new_seq]);
n_in = accumarray(seq_id, 1);
good = find(n_in == n_elem);
if numel(good) < max(seq_id)
    warning('get_abcd_trials:incomplete', '%s: dropped %d of %d sequences without %d elements (sizes %s)', ...
        paths.stimulus, max(seq_id) - numel(good), max(seq_id), n_elem, mat2str(n_in(n_in ~= n_elem)'));
end
keep = ismember(seq_id, good);
E = reshape(on(keep), n_elem, [])';
O = reshape(ori(keep), n_elem, [])';
assert(all(O == O(1, :), 'all'), 'Element orientations differ between sequences');
soa = diff(E, 1, 2);
soa_median = median(soa(:));

% element offsets
pd_on = []; pd_off = [];
if isfield(T, 'pd_on') && isfield(T.pd_on, 'sglxTime'); pd_on = T.pd_on.sglxTime(:); end
if isfield(T, 'pd_off') && isfield(T.pd_off, 'sglxTime'); pd_off = sort(T.pd_off.sglxTime(:)); end
one_pulse_per_element = numel(pd_on) == numel(on) && numel(pd_off) == numel(on);
if one_pulse_per_element
    F = nan(size(E));
    for k = 1:numel(E)
        j = find(pd_off > E(k), 1);
        if ~isempty(j) && pd_off(j) - E(k) < soa_median; F(k) = pd_off(j); end
    end
    element_dur = median(F(:) - E(:), 'omitnan');
    F(isnan(F)) = E(isnan(F)) + element_dur;
else
    element_dur = soa_median;
    F = [E(:, 2:end), E(:, end) + element_dur];
end
assert(all(abs(soa(:) - soa_median) < 0.05), 'Irregular element onsets within sequences (SOA %.3f-%.3f s)', ...
    min(soa(:)), max(soa(:)));

trials.element_names    = element_names;
trials.orientation      = O(1, :);
trials.element_onset    = E;
trials.element_offset   = F;
trials.seq_onset        = E(:, 1);
trials.seq_offset       = F(:, end);
trials.element_duration = element_dur;
trials.soa              = soa;
trials.soa_median       = soa_median;
trials.isi_median       = max(soa_median - element_dur, 0);
trials.element_start    = (0:n_elem - 1) * soa_median;
trials.element_end      = trials.element_start + element_dur;
trials.seq_duration     = trials.element_end(end);
trials.iti              = trials.seq_onset(2:end) - trials.seq_offset(1:end - 1); % grey between sequences
trials.block            = cumsum([true; trials.iti > 5]);
trials.n_dropped        = max(seq_id) - numel(good);

% every element onset should coincide with a photodiode transition
pd = sort([pd_on; pd_off]);
if isempty(pd)
    trials.pd_mismatch_max = NaN;
else
    trials.pd_mismatch_max = max(abs(interp1(pd, pd, E(:), 'nearest', 'extrap') - E(:)));
end

fprintf(['%d sequences (%d blocks) of %s, orientations %s deg, element %.0f ms every %.0f ms (grey %.0f ms ' ...
    'between elements), grey between sequences %.2f s (median), photodiode mismatch %.1f ms\n'], ...
    size(E, 1), max(trials.block), strjoin(element_names, ''), mat2str(trials.orientation), 1000 * element_dur, ...
    1000 * soa_median, 1000 * trials.isi_median, median(trials.iti), 1000 * trials.pd_mismatch_max);
end
