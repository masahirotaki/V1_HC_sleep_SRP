function plot_fsb_distribution(R, Q, out_dir, min_exp)
% PLOT_FSB_DISTRIBUTION  Distributions of spontaneous FSB (raw and rate-corrected) per unit.
%
% R = all test days of one region (both protocols pooled); Q = the broad-unit table returned
% by plot_fsb_quartiles (for the quartile edges). Units of every cell type with >= min_exp
% expected Poisson burst spikes in the spontaneous periods (as in plot_fsb_quartiles).
% Rows: all units, then L23 / L4 / L56 for V1. Columns: raw spontaneous FSB, FSB / Poisson FSB
% (log axis; units with FSB = 0 are drawn in the "0" bin at the left), spontaneous rate.
% Broad (black) and narrow (red) units overlaid as counts; triangles = medians. Lines on
% the FSB/Poisson axis: dashed = broad quartile edges over all layers (figQ1/Q2), dotted =
% broad quartile edges within that layer (figQ3).
%   figQ_fsb_distribution

if nargin < 4; min_exp = 5; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
region = "V1"; if isfield(R(1).opt, 'region'); region = string(R(1).opt.region); end
tau = R(1).opt.isi_thr;

U = table();
for i = 1:numel(R)
    u = R(i).unit;
    rate = u.n_spont / R(i).spont_dur;
    pois = 1 - exp(-2 * rate * tau);
    U = [U; table(u.cell_type, u.layer, rate, u.fsb_spont, u.fsb_spont ./ pois, u.n_spont .* pois, ...
        'VariableNames', {'cell_type', 'layer', 'rate', 'fsb', 'fsb_rc', 'exp_burst_spikes'})]; %#ok<AGROW>
end
n_all = height(U);
U = U(U.exp_burst_spikes >= min_exp, :);

row_names = "all units";
row_masks = {true(height(U), 1)};
if region == "V1"
    for l = ["L23", "L4", "L56"]
        row_names(end + 1) = l; %#ok<AGROW>
        row_masks{end + 1} = U.layer == l; %#ok<AGROW>
    end
end
edges_global = prctile(Q.fsb_rc, [25 50 75]);

zero_x = 10 ^ -2.25;                       % where FSB = 0 units are drawn on the log axis
rc_edges = 10 .^ (-2.5:0.1:1.8);
fsb_edges = 0:0.01:max(0.5, ceil(max(U.fsb) * 100) / 100);
rate_edges = 10 .^ (-1:0.1:2);
types = {'broad', U.cell_type ~= 2, [0 0 0]; 'narrow', U.cell_type == 2, [0.85 0.2 0.2]};

nR = numel(row_names);
fig = figure('Visible', 'off', 'Position', [50 50 1700 330 * nR + 80]);
tl = tiledlayout(nR, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for ri = 1:nR
    for col = 1:3
        ax = nexttile; hold on;
        leg = {};
        for ti = 1:2
            m = row_masks{ri} & types{ti, 2};
            if ~any(m); continue; end
            switch col
                case 1
                    x = U.fsb(m); e = fsb_edges;
                case 2
                    x = U.fsb_rc(m); x(x == 0) = zero_x; e = rc_edges;
                case 3
                    x = U.rate(m); e = rate_edges;
            end
            histogram(x, e, 'DisplayStyle', 'stairs', 'EdgeColor', types{ti, 3}, 'LineWidth', 1.5);
            leg{end + 1} = sprintf('%s (n=%d, median %.3g)', types{ti, 1}, nnz(m), median(x)); %#ok<AGROW>
            plot(median(x), 0, 'v', 'MarkerFaceColor', types{ti, 3}, 'MarkerEdgeColor', 'none', ...
                'MarkerSize', 8, 'HandleVisibility', 'off');
        end
        switch col
            case 1
                xlabel('spontaneous FSB (burst spikes / all spikes)');
                title(sprintf('%s: raw FSB', row_names(ri)));
            case 2
                set(ax, 'XScale', 'log');
                xline(1, 'k-', 'Poisson', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'top');
                xline(edges_global, 'b--', 'HandleVisibility', 'off');
                mb = Q.layer == row_names(ri) & ~isnan(Q.quartile_layer);
                if ri > 1 && any(mb)
                    el = arrayfun(@(k) min(Q.fsb_rc(mb & Q.quartile_layer == k)), 2:4);
                    xline(el, 'm:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
                xt = [zero_x 0.01 0.1 1 10];
                set(ax, 'XTick', xt, 'XTickLabel', {'0', '0.01', '0.1', '1', '10'});
                xlim([10 ^ -2.5, 10 ^ 1.8]);
                xlabel('FSB / Poisson FSB (log; left bin = FSB 0)');
                title(sprintf('%s: rate-corrected FSB', row_names(ri)));
            case 3
                set(ax, 'XScale', 'log');
                xlabel('spontaneous rate (Hz)');
                title(sprintf('%s: spontaneous rate', row_names(ri)));
        end
        ylabel('units');
        legend(leg, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
    end
end
title(tl, {sprintf('%s: per-unit spontaneous FSB (burst = ISI <= %g ms), test days of both protocols pooled', region, 1000 * tau), ...
    sprintf(['%d of %d units with >= %g expected Poisson burst spikes. FSB/Poisson lines: blue dashed = broad quartile edges ' ...
    'over all layers, magenta dotted = broad quartile edges within the layer'], height(U), n_all, min_exp)}, 'FontSize', 10);
save_fig_formats(fig, out_dir, 'figQ_fsb_distribution');
end
