function sessions = first_novel_sessions()
% FIRST_NOVEL_SESSIONS  One recording day per animal: the first day the animal
% sees a novel stimulus sequence, with sleep before (where recorded) and after.
% Taken from Y:\Ellie\DATA\Ellie_session_stimuli_table.xlsx. recordings are in
% acquisition order (session_id, then g number).
%
% Stimuli: GAVNIK_ABCD    = 4 static gratings, 150 ms each, back to back
%          GAVNIK250_ABCD = same, 250 ms each
%          TRAIN          = 4 gratings, 150 ms each, separated by 150 ms grey
%
% Not included: M00012 (pilot, no sleep on TRAIN day), M00070 (SparseNoise only),
% M00071 and M00098 (no Sleep_Box recordings on GAVNIK days).
% Later novel days with pre/post sleep (not listed here): M00069 20251006 (TRAIN),
% M00087 20260218 (GAVNIK_ABCD), M00088 20260316 (GAVNIK250), M00096 20260629
% (GAVNIK_ABCD), M00113 20260720 (GAVNIK250).

S = {
% subject   date        probe   novel stimulus       recordings (in order)
 'M00013', '20250203', 'NPX1', 'TRAIN',            {'OP_Tuning','TRAIN','Sleep_Box_1','Sleep_Box_2'}
 'M00014', '20250317', 'NPX1', 'GAVNIK_ABCD',      {'OP_Tuning','GAVNIK_ABCD','OP_Tuning','Sleep_Box'}
 'M00069', '20250929', 'NPX2', 'GAVNIK_ABCD_1',    {'Sleep_Box_1','GAVNIK_D___','GAVNIK_A___','GAVNIK_ABCD_1','Sleep_Box_2','GAVNIK_ABCD_2'}
 'M00087', '20260211', 'NPX2', 'GAVNIK250_ABCD_1', {'Sleep_Box_1','GAVNIK250_ABCD_1','SparseNoise','Sleep_Box_2','GAVNIK250_ABCD_2'}
 'M00088', '20260309', 'NPX2', 'GAVNIK_ABCD_1',    {'Sleep_Box_1','GAVNIK_ABCD_1','Sleep_Box_2','GAVNIK_ABCD_2'}
 'M00096', '20260622', 'NPX2', 'GAVNIK250_ABCD_1', {'Sleep_Box_1','GAVNIK250_ABCD_1','Sleep_Box_2','GAVNIK250_ABCD_2'}
 'M00100', '20260602', 'NPX2', 'GAVNIK250_ABCD_1', {'Sleep_Box_1','GAVNIK250_ABCD_1','Sleep_Box_2','GAVNIK250_ABCD_2'}
 'M00113', '20260712', 'NPX2', 'GAVNIK_ABCD_1',    {'Sleep_Box_1','GAVNIK_ABCD_1','Sleep_Box_2','GAVNIK_ABCD_2'}
};

sessions = struct('subject', S(:,1), 'date', S(:,2), 'probe_type', S(:,3), ...
    'novel_stimulus', S(:,4), 'recordings', S(:,5));

for i = 1:numel(sessions)
    rec = sessions(i).recordings;
    is_sleep = startsWith(rec, 'Sleep_Box');
    first_novel = find(strcmp(rec, sessions(i).novel_stimulus), 1);
    sessions(i).pre_sleep  = rec(is_sleep & (1:numel(rec)) < first_novel);
    sessions(i).post_sleep = rec(is_sleep & (1:numel(rec)) > first_novel);
end
end
