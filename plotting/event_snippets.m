function [snips, t_rel, kept] = event_snippets(x, tvec, fs, times, win)
% EVENT_SNIPPETS  Extract x around event times.
%   snips  [nKeptEvents x nSamples], events too close to the edges are dropped
%   t_rel  time relative to the event (s); kept = logical index into times

offsets = round(win(1) * fs):round(win(2) * fs);
centre = round((times(:) - tvec(1)) * fs) + 1;
kept = centre + offsets(1) >= 1 & centre + offsets(end) <= numel(x);
snips = reshape(x(centre(kept) + offsets), [], numel(offsets));
t_rel = offsets / fs;
end
