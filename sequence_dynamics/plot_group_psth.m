function plot_group_psth(R, out_dir, tag, fr_thr, sig_s)
% PLOT_GROUP_PSTH  All-spike PSTHs of ABCD vs A_CD vs E_CD per cell group (one protocol, sessions pooled).
%
% Columns: all units | broad, baseline < fr_thr Hz | broad, baseline >= fr_thr Hz | narrow
% (cell_type 2); baseline = mean rate from -0.6 to -0.2 s before A, averaged over the
% 3 conditions (as in plot_burst_groups). Units of all layers are pooled.
% Rows: rate (Hz / neuron, mean over neurons) and z-score (per unit, from R.unit.zmu /
% zsd, the mean / SD of single-trial 10 ms rates over the analysis window of all
% conditions). Shading = SE across neurons. sig_s = Gaussian smoothing SD in s (0 = none).
%   figP  <tag>_figP_psth_all_spikes

if nargin < 4; fr_thr = 5; end
if nargin < 5; sig_s = 0.010; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
t = R(1).t(:)';
dur = R(1).element_duration;
conds = R(1).conds;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
region = "V1"; if isfield(R(1).opt, 'region'); region = string(R(1).opt.region); end

U = pool_units(R);
C = arrayfun(@(r) double(r.psth.all), R, 'UniformOutput', false);
P = cat(1, C{:});
base = t >= -0.6 & t < -0.2;
base_fr = mean(mean(P(:, base, :), 2), 3);
P = smooth_t(P, sig_s / R(1).opt.bin);
Z = (P - U.zmu) ./ U.zsd;
active = U.zmu > 0;

broad = U.cell_type ~= 2;
groups = {true(size(broad)), broad & base_fr < fr_thr, broad & base_fr >= fr_thr, ~broad};
group_names = {'all units', sprintf('broad, base < %g Hz', fr_thr), sprintf('broad, base >= %g Hz', fr_thr), 'narrow'};
if sig_s > 0; sm_txt = sprintf('%g ms Gaussian smoothing', 1000 * sig_s); else; sm_txt = 'no smoothing'; end
mice = strjoin(arrayfun(@(r) string(r.subject(end-2:end)), R), ' ');
sess_txt = sprintf('%s units, %d ms protocol, %d mice: %s; %g ms bins, %s; shading = SE across neurons', ...
    region, R(1).protocol, numel(R), mice, 1000 * R(1).opt.bin, sm_txt);

fig = figure('Visible', 'off', 'Position', [50 50 1800 700]);
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {P, Z};
row_lab = {'rate (Hz / neuron)', 'z-score'};
for ri = 1:2
    for g = 1:4
        ax = nexttile; hold on;
        m = groups{g} & active;
        h = gobjects(1, 3);
        for c = 1:3
            h(c) = mean_se(t, rows{ri}(m, :, c), cond_col(c, :));
        end
        finish_axes(ax, t, dur, ri == 1);
        title(sprintf('%s (n=%d)', group_names{g}, nnz(m)));
        if g == 1; ylabel(row_lab{ri}); end
        if ri == 2; xlabel('Time from A/E onset (s)'); end
        if ri == 1 && g == 1; legend(h, conds, 'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none'); end
    end
end
title(tl, {'All-spike PSTHs: ABCD vs A\_CD vs E\_CD', sess_txt}, 'FontSize', 10);
save_fig_formats(fig, out_dir, [tag '_figP_psth_all_spikes']);
end


function h = mean_se(t, X, col)
mu = mean(X, 1, 'omitnan');
se = std(X, 0, 1, 'omitnan') / sqrt(max(size(X, 1), 1));
fill([t, fliplr(t)], [mu - se, fliplr(mu + se)], col, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
h = plot(t, mu, 'Color', col, 'LineWidth', 1);
end


function finish_axes(ax, t, dur, floor0)
ln = [findobj(ax, 'Type', 'line'); findobj(ax, 'Type', 'patch')];
y = cell2mat(arrayfun(@(l) l.YData(:), ln, 'UniformOutput', false));
y = y(isfinite(y));
if isempty(y); y = [0 1]; end
pad = 0.05 * max(range(y), eps);
yl = [min(y) - pad, max(y) + pad];
if floor0; yl(1) = max(yl(1), 0); end
xlim(ax, [t(1) t(end)]);
ylim(ax, yl);
p1 = patch(ax, [t(1) 0 0 t(1)], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
p2 = patch(ax, [4 * dur t(end) t(end) 4 * dur], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack([p1 p2], 'bottom');
xline(ax, 0:dur:4 * dur, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
end


function X = smooth_t(X, sig)
% Gaussian smoothing along dimension 2, renormalised at the edges (no zero-padding dip).
if sig <= 0; return; end
x = -ceil(3 * sig):ceil(3 * sig);
g = exp(-x .^ 2 / (2 * sig ^ 2));
X = convn(X, g, 'same') ./ conv(ones(1, size(X, 2)), g, 'same');
end


function U = pool_units(R)
f = fieldnames(R(1).unit);
for i = 1:numel(f)
    C = arrayfun(@(r) r.unit.(f{i}), R, 'UniformOutput', false);
    U.(f{i}) = cat(1, C{:});
end
end
