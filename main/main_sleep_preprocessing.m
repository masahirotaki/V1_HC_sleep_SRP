%% Sleep preprocessing: sleep staging, ripples (CA1), spindles and UP/DOWN states (V1 L5)
% Reads Ellie's data from Y:\Ellie\DATA and writes into each recording's analysis
% folder, Y:\Ellie\DATA\SUBJECTS\<subject>\analysis\<date>\<Sleep_Box_x>\ (see
% run_sleep_detection; only that folder and its detection\ subfolder are writable).
%
% Channels: ripples on the CA1 channel (depths_from_PSD) of shank 4 or 3 (most
%           anterior); spindles and UP/DOWN on the L5 channel of the same shank.
% Mobility: video pixel-change from extracted_behaviour.mat (Ellie's threshold).

clear; close all
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
cfg = project_config();

%% One recording
subject  = 'M00069';
date     = '20250929';
stimulus = 'Sleep_Box_1';
results = run_sleep_detection(subject, date, stimulus, cfg);

%% All pre/post sleep recordings on each animal's first novel-stimulus day
sessions = first_novel_sessions();
for i = 1:numel(sessions)
    for stimulus = [sessions(i).pre_sleep sessions(i).post_sleep]
        try
            run_sleep_detection(sessions(i).subject, sessions(i).date, stimulus{1}, cfg);
        catch err
            fprintf(2, '%s %s %s failed: %s\n', sessions(i).subject, sessions(i).date, stimulus{1}, err.message);
        end
    end
end
