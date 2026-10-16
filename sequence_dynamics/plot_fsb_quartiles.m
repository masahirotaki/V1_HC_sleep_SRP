function Q = plot_fsb_quartiles(R, out_dir, sig_s, min_exp)
% PLOT_FSB_QUARTILES  Broad units split into quartiles of rate-corrected spontaneous FSB.
%
% R = all test days of one region (both protocols). Per broad unit (cell_type ~= 2):
%   rate      spontaneous rate = r.unit.n_spont / r.spont_dur (grey periods, 3 recordings)
%   fsb_pois  FSB of a Poisson train at that rate = 1 - exp(-2 * rate * isi_thr)
%             (a spike is in a burst if the ISI before or after it is <= isi_thr)
%   fsb_rc    rate-corrected FSB = r.unit.fsb_spont / fsb_pois: 1 Poisson, < 1 refractory /
%             regular, > 1 bursty
% Units with fewer than min_exp expected Poisson burst spikes (n_spont * fsb_pois) are
% left out (fsb_rc too noisy). Quartile edges are computed on the units of both
% protocols pooled, so Q1-Q4 mean the same in the 150 and 250 ms figures.
%   figQ0_quartile_definition          fsb_rc vs rate (coloured by quartile) and per-quartile
%                                      rate / FSB / BI summaries
%   <prot>ms_figQ1_quartile_psth       columns Q1-Q4; rows: all spikes (Hz / neuron), burst
%                                      events (Hz / neuron), single spikes (Hz / neuron), burst
%                                      events / pre-stimulus baseline; ABCD vs A_CD vs E_CD
%   <prot>ms_figQ2_layer_<measure>     V1 only: one figure per row measure above, rows Q1-Q4 x
%                                      columns L23 / L4 / L56 (same quartile edges)
%   figQ0_quartile_definition_within_layer, <prot>ms_figQ3_within_layer_<measure>
%                                      V1 only: as figQ0 / figQ2 with quartiles computed within
%                                      each layer (by rank; column Q.quartile_layer)
% Returns Q, a table with one row per kept unit (subject, date, protocol, cluster_id, rate,
% fsb_spont, fsb_rc, bi_spont, quartile). sig_s = Gaussian smoothing SD in s (0 = none).

if nargin < 3; sig_s = 0.010; end
if nargin < 4; min_exp = 5; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
region = "V1"; if isfield(R(1).opt, 'region'); region = string(R(1).opt.region); end
tau = R(1).opt.isi_thr;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
q_col = [0.20 0.45 0.80; 0.45 0.70 0.35; 0.95 0.60 0.15; 0.80 0.15 0.15];

%% Per-unit rate-corrected FSB and quartiles
Q = table();
for i = 1:numel(R)
    u = R(i).unit;
    rate = u.n_spont / R(i).spont_dur;
    pois = 1 - exp(-2 * rate * tau);
    n = numel(rate);
    Q = [Q; table(repmat(string(R(i).subject), n, 1), repmat(string(R(i).date), n, 1), ...
        repmat(R(i).protocol, n, 1), repmat(i, n, 1), (1:n)', u.cluster_id, u.cell_type, u.layer, rate, ...
        u.fsb_spont, pois, u.fsb_spont ./ pois, u.bi_spont, u.n_spont .* pois, ...
        'VariableNames', {'subject', 'date', 'protocol', 'day', 'row', 'cluster_id', 'cell_type', 'layer', ...
        'rate', 'fsb_spont', 'fsb_pois', 'fsb_rc', 'bi_spont', 'exp_burst_spikes'})]; %#ok<AGROW>
end
n_broad = nnz(Q.cell_type ~= 2);
Q = Q(Q.cell_type ~= 2 & Q.exp_burst_spikes >= min_exp, :);
edges = prctile(Q.fsb_rc, [0 25 50 75 100]);
Q.quartile = discretize(Q.fsb_rc, edges, 'IncludedEdge', 'right');
q_names = arrayfun(@(k) sprintf('Q%d: FSB/Poisson %.2f-%.2f', k, edges(k), edges(k + 1)), 1:4, 'UniformOutput', false);
excl_txt = sprintf('%d of %d broad units kept (>= %g expected Poisson burst spikes in spontaneous periods)', ...
    height(Q), n_broad, min_exp);
fprintf('%s FSB/Poisson quartile edges: %s; %s\n', region, mat2str(round(edges, 2)), excl_txt);

% V1: quartiles within each layer (both protocols pooled), by rank so ties at FSB = 0
% are split evenly; units of unknown layer get NaN
layers = ["L23", "L4", "L56"];
Q.quartile_layer = nan(height(Q), 1);
edges_layer = nan(numel(layers), 5);
for li = 1:numel(layers)
    m = find(Q.layer == layers(li));
    if numel(m) < 8; continue; end
    Q.quartile_layer(m) = ceil(4 * tiedrank(Q.fsb_rc(m)) / numel(m));
    edges_layer(li, :) = [arrayfun(@(k) min(Q.fsb_rc(m(Q.quartile_layer(m) == k))), 1:4), max(Q.fsb_rc(m))];
    fprintf('%s %s within-layer quartiles (lower edges, max): %s, n = %d\n', region, layers(li), ...
        mat2str(round(edges_layer(li, :), 2)), numel(m));
end
has_layers = any(~isnan(Q.quartile_layer));

%% figQ0: quartile definition
fig = figure('Visible', 'off', 'Position', [50 50 1600 450]);
tl = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
definition_tiles(Q, Q.quartile, q_col, '');
title(tl, {sprintf('%s broad units: quartiles of rate-corrected spontaneous FSB (burst = ISI <= %g ms)', region, 1000 * tau), ...
    sprintf('%s; both protocols pooled; edges %s', excl_txt, mat2str(round(edges, 2)))}, 'FontSize', 10);
save_fig_formats(fig, out_dir, 'figQ0_quartile_definition');

if has_layers
    fig = figure('Visible', 'off', 'Position', [50 50 1600 1200]);
    tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for li = 1:numel(layers)
        m = Q.layer == layers(li) & ~isnan(Q.quartile_layer);
        definition_tiles(Q(m, :), Q.quartile_layer(m), q_col, sprintf('%s (n=%d): ', layers(li), nnz(m)));
    end
    title(tl, {sprintf('%s broad units: quartiles of rate-corrected spontaneous FSB computed WITHIN each layer', region), ...
        sprintf('%s; both protocols pooled; quartiles by rank within layer', excl_txt)}, 'FontSize', 10);
    save_fig_formats(fig, out_dir, 'figQ0_quartile_definition_within_layer');
end

%% figQ1: PSTHs per quartile, per protocol
t = R(1).t(:)';
base = t >= -0.6 & t < -0.2;
if sig_s > 0; sm_txt = sprintf('%g ms Gaussian smoothing', 1000 * sig_s); else; sm_txt = 'no smoothing'; end
rows = {'all', 'rate (Hz / neuron)', 'all spikes'; 'event', 'Hz / neuron', 'burst events'; ...
    'single', 'Hz / neuron', 'single spikes'; 'event', 'rate / baseline', 'burst events / baseline'};
for p = [150 250]
    days = find([R.protocol] == p);
    if isempty(days); continue; end
    dur = R(days(1)).element_duration;
    t = R(days(1)).t(:)';
    conds = R(days(1)).conds;
    Qp = Q(Q.protocol == p, :);
    for k = {'all', 'event', 'single'}
        X = cell(height(Qp), 1);
        for j = 1:height(Qp)
            X{j} = double(R(Qp.day(j)).psth.(k{1})(Qp.row(j), :, :));
        end
        P.(k{1}) = smooth_t(cat(1, X{:}), sig_s / R(1).opt.bin);
    end
    fig = figure('Visible', 'off', 'Position', [50 50 1800 1150]);
    tl = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for ri = 1:4
        for q = 1:4
            ax = nexttile; hold on;
            m = Qp.quartile == q;
            h = gobjects(1, 3);
            for c = 1:3
                X = P.(rows{ri, 1})(m, :, c);
                if ri == 4
                    [mu, se] = jack_norm(X, base);
                    h(c) = shade(t, mu, se, cond_col(c, :));
                else
                    h(c) = shade(t, mean(X, 1), std(X, 0, 1) / sqrt(max(nnz(m), 1)), cond_col(c, :));
                end
            end
            finish_axes(ax, t, dur, ri < 4);
            if ri == 4; yline(1, 'k:', 'HandleVisibility', 'off'); end
            if ri == 1
                title(sprintf('%s (n=%d, %d mice)\nmedian rate %.1f Hz', q_names{q}, nnz(m), ...
                    numel(unique(Qp.subject(m))), median(Qp.rate(m))));
            end
            if q == 1; ylabel({rows{ri, 3}, rows{ri, 2}}); end
            if ri == 4; xlabel('Time from A/E onset (s)'); end
            if ri == 1 && q == 1; legend(h, conds, 'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none'); end
        end
    end
    mice = strjoin(unique(extractAfter(Qp.subject, 3)), ' ');
    title(tl, {sprintf('%s broad units by rate-corrected FSB quartile, %d ms protocol (mice %s)', region, p, mice), ...
        sprintf('Rows 1-3: mean over neurons, shading = SE across neurons. Row 4: summed burst events / own baseline [-0.6, -0.2] s, jackknife SE. Burst = ISI <= %g ms; %g ms bins, %s', ...
        1000 * tau, 1000 * R(1).opt.bin, sm_txt)}, 'FontSize', 10);
    save_fig_formats(fig, out_dir, sprintf('%dms_figQ1_quartile_psth', p));

    % figQ2 (global quartile edges) / figQ3 (quartiles within layer): V1 only, one figure
    % per measure, rows Q1-Q4 x columns layer
    if ~has_layers; continue; end
    file_tag = {'all_spikes_Hz', 'event_Hz', 'single_Hz', 'event_norm'};
    versions = {'quartile', 'figQ2_layer', 'global quartile edges'; ...
                'quartile_layer', 'figQ3_within_layer', 'quartiles computed within each layer'};
    for vi = 1:2
    qv = Qp.(versions{vi, 1});
    for ri = 1:4
        fig = figure('Visible', 'off', 'Position', [50 50 1700 1250]);
        tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        for q = 1:4
            for li = 1:3
                ax = nexttile; hold on;
                m = qv == q & Qp.layer == layers(li);
                h = gobjects(1, 3);
                for c = 1:3
                    X = P.(rows{ri, 1})(m, :, c);
                    if ri == 4
                        [mu, se] = jack_norm(X, base);
                    else
                        mu = mean(X, 1); se = std(X, 0, 1) / sqrt(max(nnz(m), 1));
                    end
                    if nnz(m) == 0; mu = nan(size(t)); se = mu; end
                    h(c) = shade(t, mu, se, cond_col(c, :));
                end
                finish_axes(ax, t, dur, ri < 4);
                if ri == 4; yline(1, 'k:', 'HandleVisibility', 'off'); end
                if vi == 1
                    title(sprintf('%s, Q%d (n=%d, %d mice, median rate %.1f Hz)', layers(li), q, nnz(m), ...
                        numel(unique(Qp.subject(m))), median(Qp.rate(m))));
                    if li == 1; ylabel({q_names{q}, rows{ri, 2}}); end
                else
                    title(sprintf('%s, Q%d: FSB/Poisson %.2f-%.2f (n=%d, %d mice, %.1f Hz)', layers(li), q, ...
                        edges_layer(li, q), edges_layer(li, q + 1), nnz(m), numel(unique(Qp.subject(m))), median(Qp.rate(m))));
                    if li == 1; ylabel({sprintf('Q%d', q), rows{ri, 2}}); end
                end
                if q == 4; xlabel('Time from A/E onset (s)'); end
                if q == 1 && li == 1; legend(h, conds, 'Location', 'northwest', 'Box', 'off', 'Interpreter', 'none'); end
            end
        end
        if ri == 4
            se_txt = 'summed over neurons / own baseline [-0.6, -0.2] s, jackknife SE';
        else
            se_txt = 'mean over neurons, shading = SE across neurons';
        end
        title(tl, {sprintf('%s broad units by rate-corrected FSB quartile (%s) and layer: %s, %d ms protocol (mice %s)', ...
            region, versions{vi, 3}, rows{ri, 3}, p, mice), sprintf('%s. Burst = ISI <= %g ms; %g ms bins, %s', ...
            se_txt, 1000 * tau, 1000 * R(1).opt.bin, sm_txt)}, 'FontSize', 10);
        save_fig_formats(fig, out_dir, sprintf('%dms_%s_%s', p, versions{vi, 2}, file_tag{ri}));
    end
    end
end
end


function definition_tiles(Q, quart, q_col, prefix)
% Four tiles: FSB/Poisson vs rate coloured by quartile; rate, FSB and BI per quartile.
nexttile; hold on;
for k = 1:4
    m = quart == k;
    scatter(Q.rate(m), Q.fsb_rc(m), 10, q_col(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
yline(1, 'k:');
xlabel('spontaneous rate (Hz)'); ylabel('FSB / Poisson FSB');
title(sprintf('%srho(FSB/Poisson, rate) = %.2f', prefix, corr(Q.rate, Q.fsb_rc, 'type', 'Spearman')));
legend("Q" + (1:4), 'Location', 'best', 'Box', 'off');
vars = {'rate', 'spontaneous rate (Hz)'; 'fsb_spont', 'spontaneous FSB'; 'bi_spont', 'BI (ACG 2-6 / 20-40 ms)'};
for v = 1:3
    nexttile; hold on;
    for k = 1:4
        x = Q.(vars{v, 1})(quart == k);
        if isempty(x); continue; end
        p = prctile(x, [25 50 75]);
        plot([k k], p([1 3]), '-', 'Color', q_col(k, :), 'LineWidth', 4);
        plot(k, p(2), 'ko', 'MarkerFaceColor', 'w');
    end
    xlim([0.5 4.5]); set(gca, 'XTick', 1:4, 'XTickLabel', "Q" + (1:4));
    if v == 1; set(gca, 'YScale', 'log'); end
    ylabel(vars{v, 2}); title([vars{v, 2} ', median and IQR'], 'Interpreter', 'none');
end
end


function [y, se] = jack_norm(X, base)
% Summed rows divided by their own baseline mean, with leave-one-unit-out jackknife SE.
S = sum(X, 1);
y = S / mean(S(base));
n = size(X, 1);
if ~(mean(S(base)) > 0); y = nan(size(S)); se = y; return; end
if n < 3; se = zeros(size(y)); return; end
Si = S - X;
yi = Si ./ mean(Si(:, base), 2);
se = sqrt((n - 1) / n * sum((yi - mean(yi, 1, 'omitnan')) .^ 2, 1, 'omitnan'));
end


function h = shade(t, mu, se, col)
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
