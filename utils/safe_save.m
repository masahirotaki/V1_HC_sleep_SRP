function safe_save(file, vars, cfg)
% SAFE_SAVE  save(FILE, '-struct', VARS) after checking FILE is an allowed output
% location (check_output_path) and does not overwrite someone else's file
% (check_overwrite). Each field of the struct VARS becomes a variable in the .mat file.

check_output_path(file, cfg);
check_overwrite(file, cfg);
folder = fileparts(file);
if ~exist(folder, 'dir'); mkdir(folder); end
save(file, '-struct', 'vars', '-v7.3');
fprintf('Saved %s\n', file);
end
