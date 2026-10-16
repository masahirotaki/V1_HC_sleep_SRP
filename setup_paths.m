function repo_root = setup_paths()
% SETUP_PATHS  Add this repository (and the copied VR_NPX_analysis helpers)
% to the MATLAB path. Folders are prepended so copies here take priority over
% any VR_NPX_analysis checkout that may also be on the path.

repo_root = fileparts(mfilename('fullpath'));
all_dirs = strsplit(genpath(repo_root), pathsep);
all_dirs = all_dirs(~cellfun(@isempty, all_dirs) & ~contains(all_dirs, [filesep '.git']));
addpath(all_dirs{:});

p = strsplit(path, pathsep);
if any(contains(p, 'VR_NPX_analysis') & ~startsWith(p, repo_root))
    warning(['A VR_NPX_analysis checkout is also on the MATLAB path. ' ...
        'V1_HC_sleep_SRP copies take priority, but consider removing it to avoid confusion.']);
end
end
