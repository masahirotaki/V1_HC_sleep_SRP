function plot_seqlik_trial_structure(subject, date, stimulus_name, cfg)
% PLOT_SEQLIK_TRIAL_STRUCTURE  Where the sequence windows and the grey control windows
% of run_abcd_sequence_likelihood sit in the recording.
%
%   plot_seqlik_trial_structure('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Rows: whole recording; a zoom on the first sequences; a zoom on the first block pause;
% where grey windows start relative to the previous D offset. Grey windows that overlap
% each other are drawn on separate lanes. Saved as decoding_<stimulus>_seqlik_trial_structure.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths = get_session_paths(subject, date, stimulus_name, cfg);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
S = load(fullfile(out_dir, sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
R = S.results;
tr = R.trials;
win = R.opts.window;
dur = diff(win);
gs = R.grey_start(:);
margin = R.opts.grey_margin;
n_seq = numel(tr.seq_onset);
el_col = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655];
blue = [0.165 0.471 0.839]; grey = [0.45 0.45 0.45];

% eligible grey pool (same rule as run_abcd_sequence_likelihood)
g0 = tr.seq_offset(1:end - 1) + margin;
g1 = tr.seq_onset(2:end);
ok = g1 - g0 >= dur;
g0 = g0(ok); g1 = g1(ok);
% lanes for overlapping grey windows
[gs_sorted, o] = sort(gs);
lane = zeros(n_seq, 1);
lane_end = [];
for k = 1:n_seq
    l = find(lane_end <= gs_sorted(k), 1);
    if isempty(l); lane_end(end + 1) = 0; l = numel(lane_end); end %#ok<AGROW>
    lane(o(k)) = l;
    lane_end(l) = gs_sorted(k) + dur;
end
gap_after = arrayfun(@(s) find(tr.seq_onset < s, 1, 'last'), gs);
in_pause = (tr.seq_onset(gap_after + 1) - tr.seq_offset(gap_after)) > 5;
rel = gs - tr.seq_offset(gap_after);
n_overlap = sum(any(abs(gs - gs') < dur & ~eye(n_seq), 2));

f = figure('Color', 'w', 'Position', [30 30 1700 950]);
tl = tiledlayout(f, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s %s %s - trial structure for the sequence-likelihood analysis: %d sequence windows (%g-%g s ' ...
    'after A) and %d grey control windows (%d in short gaps, %d in block pauses; %d overlap another grey window)'], ...
    subject, date, strrep(stimulus_name, '_', ' '), n_seq, win, n_seq, sum(~in_pause), sum(in_pause), n_overlap), 'FontSize', 10.5);

% 1. whole recording
ax = nexttile(tl, [1 2]); hold(ax, 'on');
draw_layer(ax, tr, gs, lane, win, dur, g0, g1, blue, grey, false, el_col);
xlim(ax, [0 tr.seq_offset(end) + 5]);
xlabel(ax, 'time in recording (s)');
title(ax, 'whole recording: blue = sequence windows (scored), grey = grey control windows (one lane per overlap level), light band = eligible grey pool', ...
    'FontWeight', 'normal', 'FontSize', 8.5);

% 2. zoom: first sequences
ax = nexttile(tl, [1 2]); hold(ax, 'on');
draw_layer(ax, tr, gs, lane, win, dur, g0, g1, blue, grey, true, el_col);
xlim(ax, [tr.seq_onset(1) - 0.5, tr.seq_onset(6) + 0.2]);
xlabel(ax, 'time in recording (s)');
title(ax, sprintf(['zoom, sequences 1-5: A B C D (coloured), scored window %g-%g s after A (blue outline), ' ...
    'grey pool from D offset + %.1f s to the next A (light band), grey windows (dark grey)'], win, margin), ...
    'FontWeight', 'normal', 'FontSize', 8.5);

% 3. zoom: first block pause
ax = nexttile(tl, [1 2]); hold(ax, 'on');
long = find(diff(tr.seq_onset) > 5, 1);
draw_layer(ax, tr, gs, lane, win, dur, g0, g1, blue, grey, true, el_col);
if ~isempty(long)
    xlim(ax, [tr.seq_onset(long) - 3, tr.seq_onset(long + 1) + 3]);
end
xlabel(ax, 'time in recording (s)');
title(ax, 'zoom, first block pause: many grey windows land in the long pause, and they overlap', 'FontWeight', 'normal', 'FontSize', 8.5);

% 4. grey window start re previous D offset
ax = nexttile(tl); hold(ax, 'on');
histogram(ax, rel(~in_pause), 0.3:0.025:1.1, 'FaceColor', grey, 'EdgeColor', 'none');
xline(ax, margin, 'k:', 'D + margin', 'FontSize', 7);
xline(ax, median(tr.iti) - dur, 'k--', 'latest start (next A - 0.8 s)', 'FontSize', 7);
xlabel(ax, 'grey window start re previous D offset (s)'); ylabel(ax, 'windows');
title(ax, sprintf('short gaps (n = %d): every window ends at or before the next A onset', sum(~in_pause)), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);

ax = nexttile(tl); hold(ax, 'on');
histogram(ax, rel(in_pause), 25, 'FaceColor', grey, 'EdgeColor', 'none');
xlabel(ax, 'grey window start re previous D offset (s)'); ylabel(ax, 'windows');
title(ax, sprintf('block pauses (n = %d windows in %d pauses)', sum(in_pause), sum(diff(tr.seq_onset) > 5)), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);

name = sprintf('decoding_%s_seqlik_trial_structure', stimulus_name);
save_vector_figures(f, {name}, out_dir, cfg);
end


function draw_layer(ax, tr, gs, lane, win, dur, g0, g1, blue, grey, detail, el_col)
n_lane = max(lane);
y_seq = [0 1];
pool_col = [0.88 0.88 0.88];
for k = 1:numel(g0)
    patch(ax, [g0(k) g1(k) g1(k) g0(k)], [-n_lane - 0.2 -n_lane - 0.2 1.2 1.2], pool_col, 'EdgeColor', 'none');
end
for k = 1:numel(tr.seq_onset)
    a = tr.seq_onset(k) + win(1);
    if detail
        for e = 1:size(tr.element_onset, 2)
            x0 = tr.element_onset(k, e); x1 = tr.element_offset(k, e);
            patch(ax, [x0 x1 x1 x0], [0.15 0.15 0.85 0.85], el_col(e, :), 'EdgeColor', 'none');
        end
        rectangle(ax, 'Position', [a 0 dur 1], 'EdgeColor', blue, 'LineWidth', 1.2);
    else
        patch(ax, [a a + dur a + dur a], [0 0 1 1], blue, 'EdgeColor', 'none');
    end
end
for k = 1:numel(gs)
    y = -lane(k);
    patch(ax, [gs(k) gs(k) + dur gs(k) + dur gs(k)], [y y y + 0.8 y + 0.8], grey, 'EdgeColor', 'w', 'LineWidth', 0.3, ...
        'FaceAlpha', 0.85);
end
ylim(ax, [-n_lane - 0.3 1.3]);
set(ax, 'YTick', [-(n_lane:-1:1) + 0.4, 0.5], 'YTickLabel', [compose('grey lane %d', n_lane:-1:1), {'sequences'}], ...
    'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
if n_lane > 6
    yt = get(ax, 'YTick'); ytl = get(ax, 'YTickLabel');
    keep = [1:2:numel(yt) - 1, numel(yt)];
    set(ax, 'YTick', yt(unique(keep)), 'YTickLabel', ytl(unique(keep)));
end
end
