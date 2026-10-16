function check_output_path(target, cfg)
% CHECK_OUTPUT_PATH  Error unless TARGET is a recording's analysis folder, its
% detection/, Decoding/ or replay/ subfolder, or a file directly inside one of them:
%   <output_root>/<subj>/analysis/<yyyymmdd>/<stimulus>[/detection|/Decoding|/replay][/<file>]
% The <stimulus> folder must already exist (Ellie's), so nothing new is created
% elsewhere on the server. Overwriting existing files is checked in check_overwrite.

full_target = char(java.io.File(char(target)).getCanonicalPath());
full_root   = char(java.io.File(char(cfg.output_root)).getCanonicalPath());
ok = startsWith(lower(full_target), lower([full_root filesep]));
if ok
    rel = full_target(numel(full_root) + 2:end);
    parts = strsplit(rel, filesep);
    ok = numel(parts) >= 4 && strcmpi(parts{2}, 'analysis') && ~isempty(regexp(parts{3}, '^\d{8}$', 'once')) ...
        && exist(fullfile(full_root, parts{1:4}), 'dir') == 7;
    if ok && numel(parts) >= 5
        % <stimulus>/<file>, <stimulus>/<sub> or <stimulus>/<sub>/<file>, <sub> = detection | Decoding
        ok = numel(parts) == 5 || (numel(parts) == 6 && any(strcmpi(parts{5}, {'detection', 'Decoding', 'replay'})));
    end
end
if ~ok
    error('check_output_path:notAnalysisFolder', ...
        'Refusing to write "%s": outputs must be in %s\\<subj>\\analysis\\<date>\\<stimulus>[\\detection|\\Decoding|\\replay]', ...
        full_target, full_root);
end
end
