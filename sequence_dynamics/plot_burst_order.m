function plot_burst_order(R, out_dir, tag, fr_thr, sig_s)
% PLOT_BURST_ORDER  Absolute rates and recording-order checks (one protocol).
%
% Rates are the mean over neurons of each neuron's trial-averaged rate (Hz) in 5 ms
% bins, Gaussian-smoothed with SD sig_s (s; default 0 = no smoothing) before
% averaging; shading = SE across neurons (SD / sqrt(n)). Groups as in
% plot_burst_groups: layer x {narrow, broad baseline >= fr_thr Hz, broad < fr_thr}.
%   figE  burst events and single spikes (Hz), ABCD vs A_CD and ABCD vs E_CD
%   figF  first vs last R.n_sub trials of each recording (burst events, Hz)
%   figG  per mouse, broad units: burst-event and single-spike rate in the baseline,
%         A-onset and pre-B windows for the first / last trials of each recording, in
%         recording order (ABCD -> A_CD -> E_CD), normalised to the mouse's mean

if nargin < 4; fr_thr = 5; end
if nargin < 5; sig_s = 0; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
t = R(1).t(:)';
dur = R(1).element_duration;
conds = R(1).conds;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
layers = ["L23", "L4", "L56"];
group_names = ["narrow", sprintf("broad, base >= %g Hz", fr_thr), sprintf("broad, base < %g Hz", fr_thr)];
n_sub = R(1).n_sub;

U = pool_units(R);
if ~any(ismember(U.layer, layers)); layers = unique(U.layer)'; end   % e.g. CA1: one "layer"
nL = numel(layers);
fig_w = max(1250, round(1700 * nL / 3));
region = "V1"; if isfield(R(1).opt, 'region'); region = string(R(1).opt.region); end
for k = {'all', 'event', 'single'}
    C = arrayfun(@(r) double(r.psth.(k{1})), R, 'UniformOutput', false);
    P.(k{1}) = cat(1, C{:});
    C = arrayfun(@(r) double(r.psth_sub.(k{1})), R, 'UniformOutput', false);
    Psub.(k{1}) = cat(1, C{:});
end
base = t >= -0.6 & t < -0.2;
base_fr = mean(mean(P.all(:, base, :), 2), 3);
sig_bins = sig_s / R(1).opt.bin;
for k = {'event', 'single'}
    P.(k{1}) = smooth_t(P.(k{1}), sig_bins);
    Psub.(k{1}) = smooth_t(Psub.(k{1}), sig_bins);
end
if sig_s > 0; sm_txt = sprintf('%g ms Gaussian smoothing', 1000 * sig_s); else; sm_txt = 'no smoothing'; end
group = zeros(size(base_fr));
group(U.cell_type == 2) = 1;
group(U.cell_type ~= 2 & base_fr >= fr_thr) = 2;
group(U.cell_type ~= 2 & base_fr < fr_thr) = 3;
mice = strjoin(arrayfun(@(r) string(r.subject(end-2:end)), R), ' ');
sess_txt = sprintf('%s units, %d ms protocol, %d mice: %s; burst = ISI <= %g ms; %g ms bins, %s; shading = SE across neurons', ...
    region, R(1).protocol, numel(R), mice, 1000 * R(1).opt.isi_thr, 1000 * R(1).opt.bin, sm_txt);
cls_names = struct('event', 'burst events', 'single', 'single spikes');

%% figE: absolute rates, ABCD vs other
for cls = {'event', 'single'}
    for other = [2 3]
        fig = figure('Visible', 'off', 'Position', [50 50 fig_w 950]);
        tl = tiledlayout(3, nL, 'TileSpacing', 'compact', 'Padding', 'compact');
        for g = 1:3
            for li = 1:nL
                ax = nexttile; hold on;
                m = group == g & U.layer == layers(li);
                h = [mean_se(t, P.(cls{1})(m, :, 1), cond_col(1, :)), mean_se(t, P.(cls{1})(m, :, other), cond_col(other, :))];
                finish_axes(ax, t, dur);
                title(sprintf('%s, %s (n=%d)', layers(li), group_names(g), nnz(m)));
                if li == 1; ylabel(sprintf('%s rate (Hz / neuron)', cls_names.(cls{1}))); end
                if g == 3; xlabel('Time from A/E onset (s)'); end
                if g == 1 && li == 1; legend(h, conds([1 other]), 'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none'); end
            end
        end
        title(tl, {sprintf('%s, %s vs %s', cls_names.(cls{1}), conds{1}, conds{other}), sess_txt}, 'FontSize', 10, 'Interpreter', 'none');
        save_fig_formats(fig, out_dir, sprintf('%s_figE_%s_Hz_%s_vs_%s', tag, cls{1}, conds{1}, conds{other}));
    end
end

%% figF: first vs last trials of each recording (burst events)
for other = [2 3]
    fig = figure('Visible', 'off', 'Position', [50 50 fig_w 950]);
    tl = tiledlayout(3, nL, 'TileSpacing', 'compact', 'Padding', 'compact');
    for g = 1:3
        for li = 1:nL
            ax = nexttile; hold on;
            m = group == g & U.layer == layers(li);
            h = gobjects(1, 4); k = 0;
            for c = [1 other]
                for s = 1:2
                    k = k + 1;
                    col = cond_col(c, :);
                    if s == 2; col = 0.55 * col + 0.45; end   % last trials: lighter
                    h(k) = mean_se(t, Psub.event(m, :, c, s), col);
                    if s == 2; h(k).LineStyle = '--'; end
                end
            end
            finish_axes(ax, t, dur);
            title(sprintf('%s, %s (n=%d)', layers(li), group_names(g), nnz(m)));
            if li == 1; ylabel('burst events (Hz / neuron)'); end
            if g == 3; xlabel('Time from A/E onset (s)'); end
            if g == 1 && li == 1
                legend(h, {sprintf('%s first %d', conds{1}, n_sub), sprintf('%s last %d', conds{1}, n_sub), ...
                    sprintf('%s first %d', conds{other}, n_sub), sprintf('%s last %d', conds{other}, n_sub)}, ...
                    'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none');
            end
        end
    end
    title(tl, {sprintf('Burst events: first vs last %d trials of each recording (recorded %s then %s)', ...
        n_sub, conds{1}, conds{other}), sess_txt}, 'FontSize', 10, 'Interpreter', 'none');
    save_fig_formats(fig, out_dir, sprintf('%s_figF_first_last_%s_vs_%s', tag, conds{1}, conds{other}));
end

%% figG: per-mouse windows in recording order
wins = [1 2 4];   % R.trial_win rows: baseline, A onset, pre-B
win_txt = arrayfun(@(k) sprintf('[%.2f, %.2f] s', R(1).trial_win(k, :)), wins, 'UniformOutput', false);
win_txt = strcat({'baseline ', 'A onset ', 'pre-B '}, win_txt);
pos_names = {};
for c = 1:3; pos_names = [pos_names, {sprintf('%s 1-%d', conds{c}, n_sub), sprintf('%s last %d', conds{c}, n_sub)}]; end %#ok<AGROW>
fig = figure('Visible', 'off', 'Position', [50 50 1700 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
mouse_col = lines(numel(R));
for ci = 1:2
    cls = {'event', 'single'};
    for wi = 1:3
        nexttile; hold on;
        V = nan(numel(R), 6);
        for s = 1:numel(R)
            m = R(s).unit.cell_type ~= 2;
            wdur = diff(R(s).trial_win(wins(wi), :));
            for c = 1:3
                X = R(s).trial_counts.(cls{ci}){c}(m, :, wins(wi));
                nT = size(X, 2);
                V(s, 2 * c - 1) = mean(X(:, 1:n_sub), 'all') / wdur;
                V(s, 2 * c) = mean(X(:, nT - n_sub + 1:nT), 'all') / wdur;
            end
        end
        Vn = V ./ mean(V, 2);
        for s = 1:numel(R)
            plot(1:6, Vn(s, :), '-o', 'Color', mouse_col(s, :), 'MarkerSize', 4);
        end
        errorbar(1:6, mean(Vn, 1), std(Vn, 0, 1) / sqrt(numel(R)), 'k-', 'LineWidth', 2);
        xline([2.5 4.5], 'k:');
        yline(1, 'k:');
        set(gca, 'XTick', 1:6, 'XTickLabel', pos_names, 'TickLabelInterpreter', 'none');
        xtickangle(30);
        xlim([0.5 6.5]);
        title(sprintf('%s, %s', cls_names.(cls{ci}), win_txt{wi}));
        if wi == 1; ylabel('rate / mouse mean (broad units)'); end
        if ci == 1 && wi == 1
            legend([arrayfun(@(r) string(r.subject), R), "mean +- SE"], 'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
    end
end
title(tl, {sprintf('Recording-order check: first / last %d trials of each recording, in recording order', n_sub), ...
    sess_txt}, 'FontSize', 10);
save_fig_formats(fig, out_dir, sprintf('%s_figG_order_per_mouse', tag));
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
% Gaussian smoothing along dimension 2 (sig in bins; <= 0 leaves X unchanged).
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
