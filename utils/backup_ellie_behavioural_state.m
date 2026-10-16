function backup_ellie_behavioural_state(folder, cfg)
% BACKUP_ELLIE_BEHAVIOURAL_STATE  Keep Ellie's original behavioural_state.mat as
% behavioural_state_Ellie.mat before ours replaces it. Ours is recognised by its
% .info field. The backup name is not in cfg.output.own_files, so it can never be
% overwritten afterwards.

original = fullfile(folder, 'behavioural_state.mat');
backup = fullfile(folder, 'behavioural_state_Ellie.mat');
if ~isfile(original) || isfile(backup); return; end
S = load(original);
if isfield(S, 'behavioural_state') && isfield(S.behavioural_state, 'info'); return; end % already ours
check_output_path(backup, cfg);
copyfile(original, backup);
a = dir(original);
b = dir(backup);
assert(numel(b) == 1 && b.bytes == a.bytes, 'Backup of %s failed', original);
fprintf('Backed up Ellie''s behavioural_state.mat -> %s\n', backup);
end
