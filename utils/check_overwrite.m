function check_overwrite(target_file, cfg)
% CHECK_OVERWRITE  Error if TARGET_FILE exists and is not one of this pipeline's outputs.
% Allowed names are cfg.output.own_files.<subfolder> for files in a known subfolder
% (e.g. detection/), otherwise cfg.output.own_files.stimulus.

if exist(target_file, 'file') ~= 2; return; end
[folder, name, ext] = fileparts(target_file);
[~, parent] = fileparts(folder);
key = lower(parent);
if isfield(cfg.output.own_files, key)
    patterns = cfg.output.own_files.(key);
else
    patterns = cfg.output.own_files.stimulus;
end
regex = regexptranslate('wildcard', patterns);
if ~any(cellfun(@(r) ~isempty(regexpi([name ext], ['^' r '$'], 'once')), regex))
    error('check_overwrite:wouldOverwrite', ...
        'Refusing to overwrite existing file "%s" (not one of this pipeline''s outputs).', target_file);
end
end
