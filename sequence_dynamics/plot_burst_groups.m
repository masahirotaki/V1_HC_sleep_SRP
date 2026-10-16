function plot_burst_groups(R, out_dir, tag, fr_thr, sig_s)
% PLOT_BURST_GROUPS  Burst time courses per layer x cell group (one protocol, sessions pooled).
%
% Groups: narrow (cell_type 2), broad with baseline rate >= fr_thr Hz, broad < fr_thr;
% baseline = mean rate from -0.6 to -0.2 s before A, averaged over the 3 conditions.
% Units are weighted equally: each unit's counts are divided by its mean rate over
% the analysis window, so high-rate units do not dominate the pooled ratios.
% Shading = jackknife SE across units (leave one unit out, recompute the curve).
% sig_s = Gaussian smoothing SD in s (default 0.010; 0 = no smoothing, 5 ms bins).
%   figA  FSB(t) (burst spikes / all spikes), ABCD vs A_CD vs E_CD
%   figB  burst events and single spikes / pre-stimulus baseline, ABCD vs A_CD
%   figC  same, ABCD vs E_CD
%   figD  diagnostics for between-recording differences: per-session FSB in the
%         baseline / A / pre-B windows by condition, and FSB across trials

if nargin < 4; fr_thr = 10; end
if nargin < 5; sig_s = 0.010; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
t = R(1).t(:)';
dur = R(1).element_duration;
bin = R(1).opt.bin;
conds = R(1).conds;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
layers = ["L23", "L4", "L56"];
group_names = ["narrow", sprintf("broad, base >= %g Hz", fr_thr), sprintf("broad, base < %g Hz", fr_thr)];
thr_ms = 1000 * R(1).opt.isi_thr;

U = pool_units(R);
if ~any(ismember(U.layer, layers)); layers = unique(U.layer)'; end   % e.g. CA1: one "layer"
nL = numel(layers);
fig_w = max(1250, round(1700 * nL / 3));
region = "V1"; if isfield(R(1).opt, 'region'); region = string(R(1).opt.region); end
for k = {'all', 'burst', 'event', 'single'}
    C = arrayfun(@(r) double(r.psth.(k{1})), R, 'UniformOutput', false);
    P.(k{1}) = cat(1, C{:});
end
base = t >= -0.6 & t < -0.2;
base_fr = mean(mean(P.all(:, base, :), 2), 3);
mean_fr = mean(mean(P.all, 2), 3);
w = 1 ./ max(mean_fr, 0.1);
w(mean_fr < 0.2) = 0;
group = zeros(size(base_fr));
group(U.cell_type == 2) = 1;
group(U.cell_type ~= 2 & base_fr >= fr_thr) = 2;
group(U.cell_type ~= 2 & base_fr < fr_thr) = 3;
sm = @(x) smooth_t(x, sig_s / bin);
wrows = @(X, m, c) X(m, :, c) .* w(m);
if sig_s > 0; sm_txt = sprintf('%g ms Gaussian smoothing', 1000 * sig_s); else; sm_txt = 'no smoothing'; end
sess_txt = sprintf('%s units, %d ms protocol, %d mice: %s; units weighted equally; burst = ISI <= %g ms; %g ms bins, %s; shading = jackknife SE across units', ...
    region, R(1).protocol, numel(R), strjoin(arrayfun(@(r) string(r.subject(end-2:end)), R), ' '), thr_ms, 1000 * bin, sm_txt);

%% figA: FSB(t)
fig = figure('Visible', 'off', 'Position', [50 50 fig_w 950]);
tl = tiledlayout(3, nL, 'TileSpacing', 'compact', 'Padding', 'compact');
for g = 1:3
    for li = 1:nL
        ax = nexttile; hold on;
        m = group == g & U.layer == layers(li) & w > 0;
        h = gobjects(1, 3);
        for c = 1:3
            [y, se] = jack_ratio(sm(wrows(P.burst, m, c)), sm(wrows(P.all, m, c)));
            h(c) = shade_line(t, y, se, cond_col(c, :), '-');
        end
        finish_axes(ax, t, dur);
        title(sprintf('%s, %s (n=%d)', layers(li), group_names(g), nnz(m)));
        if li == 1; ylabel('FSB'); end
        if g == 3; xlabel('Time from A/E onset (s)'); end
        if g == 1 && li == 1; legend(h, strrep(conds, '_', '\_'), 'Location', 'northwest', 'Box', 'off'); end
    end
end
title(tl, {'Fraction of spikes in bursts over the trial', sess_txt}, 'FontSize', 10);
save_fig_formats(fig, out_dir, [tag '_figA_fsb_groups']);

%% figB / figC: burst events and singles relative to baseline
for other = [2 3]
    fig = figure('Visible', 'off', 'Position', [50 50 fig_w 950]);
    tl = tiledlayout(3, nL, 'TileSpacing', 'compact', 'Padding', 'compact');
    for g = 1:3
        for li = 1:nL
            ax = nexttile; hold on;
            m = group == g & U.layer == layers(li) & w > 0;
            h = gobjects(1, 4);
            for c = [1 other]
                ls = '-'; if c ~= 1; ls = '--'; end
                [y, se] = jack_norm(sm(wrows(P.event, m, c)), base);
                h(1 + 2 * (c ~= 1)) = shade_line(t, y, se, [0.80 0.10 0.10], ls);
                [y, se] = jack_norm(sm(wrows(P.single, m, c)), base);
                h(2 + 2 * (c ~= 1)) = shade_line(t, y, se, [0.10 0.10 0.10], ls);
            end
            finish_axes(ax, t, dur);
            yline(1, 'k:', 'HandleVisibility', 'off');
            title(sprintf('%s, %s (n=%d)', layers(li), group_names(g), nnz(m)));
            if li == 1; ylabel('rate / pre-stimulus baseline'); end
            if g == 3; xlabel('Time from A/E onset (s)'); end
            if g == 1 && li == 1
                legend(h, {['burst events ' conds{1}], ['single spikes ' conds{1}], ...
                    ['burst events ' conds{other}], ['single spikes ' conds{other}]}, ...
                    'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
            end
        end
    end
    title(tl, {sprintf('Burst events (red) and single spikes (black): %s solid vs %s dashed', ...
        conds{1}, conds{other}), sess_txt}, 'FontSize', 10, 'Interpreter', 'none');
    save_fig_formats(fig, out_dir, sprintf('%s_fig%s_events_vs_singles_%s_vs_%s', tag, char('A' + other - 1), conds{1}, conds{other}));
end

%% figD: diagnostics (broad units)
tw_sel = [1 3 4];   % R.trial_win rows: baseline, A element, pre-B
win_names = arrayfun(@(k) sprintf('[%.2f, %.2f] s', R(1).trial_win(k, :)), tw_sel, 'UniformOutput', false);
win_names = strcat({'baseline ', 'A element ', 'pre-B '}, win_names);
fig = figure('Visible', 'off', 'Position', [50 50 1600 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for wi = 1:3
    nexttile; hold on;
    leg = strings(1, numel(R));
    for s = 1:numel(R)
        m = R(s).unit.cell_type ~= 2;
        f = zeros(1, 3);
        for c = 1:3
            b = R(s).trial_counts.burst{c}(m, :, tw_sel(wi)); a = R(s).trial_counts.all{c}(m, :, tw_sel(wi));
            f(c) = sum(b(:)) / max(sum(a(:)), 1);
        end
        plot(1:3, f, '-o', 'LineWidth', 1.2);
        [~, ord] = sort(double(R(s).rec_order(3, :)));
        leg(s) = sprintf('%s (order %s)', R(s).subject, strjoin(string(ord), ''));
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', strrep(conds, '_', '\_'));
    xlim([0.7 3.3]); ylabel('FSB, broad units'); title(win_names{wi});
    if wi == 1; legend(leg, 'Location', 'best', 'Box', 'off', 'FontSize', 7); end
end
n_bin = 20;
for wi = 1:3
    nexttile; hold on;
    for c = 1:3
        bsum = 0; asum = 0;
        nT_min = min(arrayfun(@(r) r.n_trials(c), R));
        for s = 1:numel(R)
            m = R(s).unit.cell_type ~= 2;
            bsum = bsum + squeeze(sum(R(s).trial_counts.burst{c}(m, 1:nT_min, tw_sel(wi)), 1));
            asum = asum + squeeze(sum(R(s).trial_counts.all{c}(m, 1:nT_min, tw_sel(wi)), 1));
        end
        nT = numel(bsum);
        edges = 0:n_bin:nT;
        k = discretize(1:nT, edges + 0.5);
        f = accumarray(k(:), bsum(:)) ./ accumarray(k(:), asum(:));
        plot(edges(1:end - 1) + n_bin / 2, f, '-o', 'Color', cond_col(c, :), 'LineWidth', 1.2);
    end
    xline(50.5:50:150.5, 'k:');
    xlabel(sprintf('Trial (%d-trial bins; dotted = block pauses)', n_bin)); ylabel('FSB, broad units');
    title(win_names{wi});
    if wi == 1; legend(strrep(conds, '_', '\_'), 'Location', 'best', 'Box', 'off'); end
end
title(tl, {'Between-recording diagnostics. Top: per session (order = recording order of ABCD, A\_CD, E\_CD); bottom: pooled across sessions over trials', ...
    sess_txt}, 'FontSize', 10);
save_fig_formats(fig, out_dir, [tag '_figD_recording_diagnostics']);
end


function [y, se] = jack_ratio(num, den)
% Pooled ratio sum(num)/sum(den) over unit rows, with leave-one-unit-out jackknife SE.
y = sum(num, 1) ./ sum(den, 1);
n = size(num, 1);
if n < 3; se = zeros(size(y)); return; end
yi = (sum(num, 1) - num) ./ (sum(den, 1) - den);
se = sqrt((n - 1) / n * sum((yi - mean(yi, 1, 'omitnan')) .^ 2, 1, 'omitnan'));
end


function [y, se] = jack_norm(X, base)
% Pooled sum over unit rows divided by its own mean in the baseline bins, with jackknife SE.
S = sum(X, 1);
y = S / mean(S(base));
n = size(X, 1);
if n < 3; se = zeros(size(y)); return; end
Si = S - X;
yi = Si ./ mean(Si(:, base), 2);
se = sqrt((n - 1) / n * sum((yi - mean(yi, 1)) .^ 2, 1));
end


function h = shade_line(t, y, se, col, ls)
fill([t, fliplr(t)], [y - se, fliplr(y + se)], col, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
h = plot(t, y, ls, 'Color', col, 'LineWidth', 1.2);
end


function finish_axes(ax, t, dur)
ln = [findobj(ax, 'Type', 'line'); findobj(ax, 'Type', 'patch')];
y = cell2mat(arrayfun(@(l) l.YData(:), ln, 'UniformOutput', false));
y = y(isfinite(y));
if isempty(y); y = [0 1]; end
pad = 0.05 * max(range(y), eps);
yl = [min(y) - pad, max(y) + pad];
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


function U = pool_units(R)
f = fieldnames(R(1).unit);
for i = 1:numel(f)
    C = arrayfun(@(r) r.unit.(f{i}), R, 'UniformOutput', false);
    U.(f{i}) = cat(1, C{:});
end
end
