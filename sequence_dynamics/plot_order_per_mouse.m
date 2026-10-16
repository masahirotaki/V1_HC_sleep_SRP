function plot_order_per_mouse(r, out_dir, sig_s)
% PLOT_ORDER_PER_MOUSE  Recording-order check for one test day (one mouse, one protocol).
%
% Rows 1-2: burst-event and single-spike rate (Hz / neuron) in the R.trial_win windows
%   (baseline, A onset, A element, pre-B) for the first / last r.n_sub trials of each
%   recording, in recording order; broad (black) and narrow (grey) units, mean +- SE
%   across neurons.
% Row 3: burst-event PSTH of broad units (Hz / neuron), first vs last r.n_sub trials of
%   each recording, Gaussian-smoothed with SD sig_s (s; 0 = none).

if nargin < 3; sig_s = 0; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
t = r.t(:)';
dur = r.element_duration;
conds = r.conds;
n_sub = r.n_sub;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
groups = {r.unit.cell_type ~= 2, r.unit.cell_type == 2};
group_names = {'broad', 'narrow'};
group_col = [0 0 0; 0.6 0.6 0.6];
win_names = {'baseline', 'A onset', 'A element', 'pre-B'};
cls = {'event', 'single'};
cls_names = {'burst events', 'single spikes'};
pos_names = {};
for c = 1:3; pos_names = [pos_names, {sprintf('%s 1-%d', conds{c}, n_sub), sprintf('%s last %d', conds{c}, n_sub)}]; end %#ok<AGROW>

fig = figure('Visible', 'off', 'Position', [50 50 1800 1100]);
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for ci = 1:2
    for wi = 1:size(r.trial_win, 1)
        nexttile; hold on;
        wdur = diff(r.trial_win(wi, :));
        h = gobjects(1, 2);
        for g = 1:2
            m = groups{g};
            if ~any(m); continue; end
            V = nan(nnz(m), 6);
            for c = 1:3
                X = double(r.trial_counts.(cls{ci}){c}(m, :, wi));
                nT = size(X, 2);
                V(:, 2 * c - 1) = mean(X(:, 1:n_sub), 2) / wdur;
                V(:, 2 * c) = mean(X(:, nT - n_sub + 1:nT), 2) / wdur;
            end
            h(g) = errorbar(1:6, mean(V, 1), std(V, 0, 1) / sqrt(size(V, 1)), '-o', ...
                'Color', group_col(g, :), 'LineWidth', 1.5, 'MarkerFaceColor', group_col(g, :));
        end
        xline([2.5 4.5], 'k:');
        set(gca, 'XTick', 1:6, 'XTickLabel', pos_names, 'TickLabelInterpreter', 'none');
        xtickangle(30); xlim([0.5 6.5]);
        title(sprintf('%s, %s [%.2f, %.2f] s', cls_names{ci}, win_names{wi}, r.trial_win(wi, :)));
        if wi == 1; ylabel('Hz / neuron'); end
        if ci == 1 && wi == 1
            ok = isgraphics(h);
            leg = cellfun(@(nm, m) sprintf('%s (n=%d)', nm, nnz(m)), group_names(ok), groups(ok), 'UniformOutput', false);
            legend(h(ok), leg, 'Location', 'best', 'Box', 'off');
        end
    end
end
m = groups{1};
for c = 1:3
    ax = nexttile; hold on;
    for s = 1:2
        col = cond_col(c, :);
        if s == 2; col = 0.55 * col + 0.45; end
        X = smooth_t(double(r.psth_sub.event(m, :, c, s)), sig_s / r.opt.bin);
        hh = mean_se(t, X, col);
        if s == 2; hh.LineStyle = '--'; end
    end
    finish_axes(ax, t, dur);
    legend({sprintf('first %d', n_sub), sprintf('last %d', n_sub)}, 'Location', 'northeast', 'Box', 'off');
    title(sprintf('%s: broad-unit burst events', conds{c}), 'Interpreter', 'none');
    xlabel('Time from A/E onset (s)');
    if c == 1; ylabel('Hz / neuron'); end
end
ax = nexttile; hold on;
for c = 1:3
    mean_se(t, smooth_t(double(r.psth.event(m, :, c)), sig_s / r.opt.bin), cond_col(c, :));
end
finish_axes(ax, t, dur);
legend(conds, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
title('all trials, broad-unit burst events');
xlabel('Time from A/E onset (s)');

if sig_s > 0; sm_txt = sprintf('%g ms Gaussian smoothing', 1000 * sig_s); else; sm_txt = 'no smoothing'; end
order_txt = strjoin(r.rec_order(1, :) + " (" + r.rec_order(2, :) + "/g" + r.rec_order(3, :) + ")", ' -> ');
region = "V1"; if isfield(r.opt, 'region'); region = string(r.opt.region); end
title(tl, {sprintf('%s %s, %s units, %d ms protocol (week %d): recording-order check', r.subject, r.date, region, r.protocol, r.week), ...
    sprintf('Recorded: %s. Burst = ISI <= %g ms; PSTHs %g ms bins, %s; mean +- SE across neurons', ...
    order_txt, 1000 * r.opt.isi_thr, 1000 * r.opt.bin, sm_txt)}, 'FontSize', 10, 'Interpreter', 'none');
save_fig_formats(fig, out_dir, sprintf('%s_%s_%dms_order_check', r.subject, r.date, r.protocol));
end


function h = mean_se(t, X, col)
mu = mean(X, 1, 'omitnan');
se = std(X, 0, 1, 'omitnan') / sqrt(max(size(X, 1), 1));
fill([t, fliplr(t)], [mu - se, fliplr(mu + se)], col, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
h = plot(t, mu, 'Color', col, 'LineWidth', 1);
end


function finish_axes(ax, t, dur)
ln = [findobj(ax, 'Type', 'line'); findobj(ax, 'Type', 'patch')];
y = cell2mat(arrayfun(@(l) l.YData(:), ln, 'UniformOutput', false));
y = y(isfinite(y));
if isempty(y); y = [0 1]; end
pad = 0.05 * max(range(y), eps);
yl = [max(min(y) - pad, 0), max(y) + pad];
xlim(ax, [t(1) t(end)]);
ylim(ax, yl);
p1 = patch(ax, [t(1) 0 0 t(1)], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
p2 = patch(ax, [4 * dur t(end) t(end) 4 * dur], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack([p1 p2], 'bottom');
xline(ax, 0:dur:4 * dur, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
end


function X = smooth_t(X, sig)
if sig <= 0; return; end
x = -ceil(3 * sig):ceil(3 * sig);
g = exp(-x .^ 2 / (2 * sig ^ 2));
X = convn(X, g / sum(g), 'same');
end
