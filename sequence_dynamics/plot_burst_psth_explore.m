function plot_burst_psth_explore(R, out_dir, tag, compact)
% PLOT_BURST_PSTH_EXPLORE  Exploratory figures from compute_session_burst_psth (one protocol).
%
% R may hold several sessions of the same protocol; their units are pooled.
% Broad = cell_type ~= 2 (peak-to-valley > 0.425 ms), narrow = cell_type == 2.
% Burst index (BI) = unit.bi_spont (spontaneous ACG 2-6 ms vs 20-40 ms, log2); terciles
% over broad units with >= 20 spontaneous ACG pairs at 20-40 ms.
%   fig1  z-scored PSTH per condition, layer x BI tercile (+ narrow units)
%   fig2  time course of the fraction of spikes in bursts; burst events vs single spikes
%   fig3  unit heatmaps sorted by layer then BI   (compact = true: fig1 and fig3 only)
%   fig4  burstiness characterisation (depth, rate, FSB, LV, ACG)
%   fig5  per-unit epoch metrics vs BI
%   fig6  k-means clusters of ABCD response shape, with layer / BI composition

if nargin < 4; compact = false; end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
assert(isscalar(unique([R.protocol])), 'Pool one protocol at a time');
t = R(1).t(:)';
dur = R(1).element_duration;
bin = R(1).opt.bin;
conds = R(1).conds;
cond_col = [0.25 0.25 0.85; 0.90 0.25 0.65; 0.15 0.70 0.60];
layers = ["L23", "L4", "L56"];
layer_col = [0.85 0.40 0.10; 0.20 0.60 0.20; 0.30 0.30 0.80];
bi_names = {'low', 'mid', 'high'};
bi_col = [0.55 0.55 0.55; 0.95 0.65 0.20; 0.80 0.10 0.10];

U = pool_units(R);
for k = {'all', 'burst', 'event', 'single'}
    C = arrayfun(@(r) double(r.psth.(k{1})), R, 'UniformOutput', false);
    P.(k{1}) = cat(1, C{:});
end
sm = @(X) smooth_t(X, 0.010 / bin);
Z = (sm(P.all) - U.zmu) ./ U.zsd;

broad = U.cell_type ~= 2;
narrow = U.cell_type == 2;
ok_bi = broad & U.acg_long >= 20;
q = quantile(U.bi_spont(ok_bi), [1/3 2/3]);
bicls = zeros(size(broad));
bicls(ok_bi & U.bi_spont <= q(1)) = 1;
bicls(ok_bi & U.bi_spont > q(1) & U.bi_spont <= q(2)) = 2;
bicls(ok_bi & U.bi_spont > q(2)) = 3;
sess_txt = sprintf('%d ms protocol, %d session(s): %s', R(1).protocol, numel(R), ...
    strjoin(arrayfun(@(r) [r.subject ' ' r.date], R, 'UniformOutput', false), ', '));

%% fig1: PSTH by layer x burstiness
fig = figure('Visible', 'off', 'Position', [50 50 1600 950]);
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for li = 1:3
    for col = 1:4
        ax = nexttile;
        if col <= 3
            m = U.layer == layers(li) & bicls == col;
            ttl = sprintf('%s broad, BI %s (n=%d)', layers(li), bi_names{col}, nnz(m));
        else
            m = U.layer == layers(li) & narrow;
            ttl = sprintf('%s narrow (n=%d)', layers(li), nnz(m));
        end
        h = plot_conditions(t, Z(m, :, :), dur, cond_col);
        title(ttl);
        if li == 3; xlabel('Time from A/E onset (s)'); end
        if col == 1; ylabel('FR (z)'); end
        if li == 1 && col == 1
            legend(ax, h, strrep(conds, '_', '\_'), 'Location', 'northwest', 'Box', 'off');
        end
    end
end
title(tl, {sprintf('V1 PSTH by layer and spontaneous burst index (BI terciles %.2f / %.2f)', q), sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig1_psth_layer_by_burstiness']);

%% fig3: heatmaps sorted by layer then BI
m = find(broad & U.layer ~= "unknown");
[~, li_idx] = ismember(U.layer(m), layers);
bi_sort = U.bi_spont(m);
bi_sort(~ok_bi(m)) = NaN;
[~, o] = sortrows([li_idx, -bi_sort]);
idx = m(o);
bounds = find(diff(li_idx(o))) + 0.5;
fig = figure('Visible', 'off', 'Position', [50 50 1700 900]);
tl = tiledlayout(1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {Z(idx, :, 1), Z(idx, :, 2), Z(idx, :, 3), Z(idx, :, 2) - Z(idx, :, 3), Z(idx, :, 1) - Z(idx, :, 2)};
names = {conds{1}, conds{2}, conds{3}, [conds{2} ' - ' conds{3}], [conds{1} ' - ' conds{2}]};
for k = 1:numel(panels)
    nexttile;
    imagesc(t, 1:numel(idx), panels{k});
    if k <= 3; clim([-1 3]); colormap(gca, parula); else; clim([-1.5 1.5]); colormap(gca, diverging_map()); end
    hold on;
    xline(0:dur:4 * dur, 'w-', 'LineWidth', 0.5);
    yline(bounds, 'k-', 'LineWidth', 1.5);
    title(strrep(names{k}, '_', '\_'));
    xlabel('Time (s)');
    if k == 1
        ylabel(sprintf('Broad V1 units: %s (top to bottom), high BI first', strjoin(layers, ' / ')));
    else
        set(gca, 'YTickLabel', []);
    end
    colorbar;
end
nexttile;
imagesc(1, 1:numel(idx), bi_sort(o), 'AlphaData', ~isnan(bi_sort(o)));
clim(quantile(U.bi_spont(ok_bi), [0.05 0.95]));
colormap(gca, hot);
yline(bounds, 'c-', 'LineWidth', 1.5);
set(gca, 'XTick', [], 'YTickLabel', [], 'Color', [0.8 0.8 0.8]);
title('BI (grey = too few spikes)');
colorbar;
title(tl, {'Unit z-scored PSTHs', sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig3_heatmaps']);

if compact; return; end

%% fig2: time course of burstiness
fig = figure('Visible', 'off', 'Position', [50 50 1500 1000]);
tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
row_names = [layers + " broad", "all narrow"];
base = t < -0.2;
for li = 1:4
    if li <= 3; m = broad & U.layer == layers(li); else; m = narrow; end
    % pooled fraction of spikes in bursts, per condition
    ax = nexttile; hold on;
    h = gobjects(1, 3);
    for c = 1:3
        fsb = sm(sum(P.burst(m, :, c), 1)) ./ sm(sum(P.all(m, :, c), 1));
        h(c) = plot(t, fsb, 'Color', cond_col(c, :), 'LineWidth', 1.2);
    end
    finish_axes(ax, t, dur);
    title(sprintf('%s (n=%d): fraction of spikes in bursts', row_names(li), nnz(m)));
    ylabel('FSB');
    if li == 1; legend(h, strrep(conds, '_', '\_'), 'Location', 'best', 'Box', 'off'); end
    % burst events vs single spikes, ABCD and E_CD
    ax = nexttile; hold on;
    styles = {'-', ':'};
    for c = [1 3]
        ev = sm(sum(P.event(m, :, c), 1)); sg = sm(sum(P.single(m, :, c), 1));
        s = styles{1 + (c == 3)};
        plot(t, ev / mean(ev(base)), s, 'Color', [0.8 0.1 0.1], 'LineWidth', 1.2);
        plot(t, sg / mean(sg(base)), s, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
    end
    finish_axes(ax, t, dur);
    yline(1, 'k:');
    title('burst events (red) vs single spikes (black); ABCD solid, E\_CD dotted');
    ylabel('rate / baseline');
    % ABCD FSB by BI tercile
    ax = nexttile; hold on;
    if li <= 3
        h = gobjects(1, 3);
        for b = 1:3
            mb = m & bicls == b;
            fsb = sm(sum(P.burst(mb, :, 1), 1)) ./ sm(sum(P.all(mb, :, 1), 1));
            h(b) = plot(t, fsb, 'Color', bi_col(b, :), 'LineWidth', 1.2);
        end
        finish_axes(ax, t, dur);
        title('ABCD fraction of spikes in bursts, by BI tercile');
        if li == 1; legend(h, strcat("BI ", bi_names), 'Location', 'best', 'Box', 'off'); end
    else
        axis off;
    end
    if li >= 3; xlabel('Time from A/E onset (s)'); end
end
title(tl, {'Burstiness over the trial (pooled spikes, 10 ms Gaussian smoothing)', sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig2_burst_timecourse']);

%% fig4: burstiness characterisation
fig = figure('Visible', 'off', 'Position', [50 50 1500 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
fr = U.fr(:, 1);
ok_n = narrow & U.acg_long >= 20;
nexttile; hold on;
for li = 1:3
    m = ok_bi & U.layer == layers(li);
    scatter(U.bi_spont(m), U.rel_L4(m), 12, layer_col(li, :), 'filled', 'MarkerFaceAlpha', 0.6);
end
scatter(U.bi_spont(ok_n), U.rel_L4(ok_n), 12, 'k', 'x');
yline([-60 60], 'k:'); xline(0, 'k:');
xlabel('BI (spontaneous ACG)'); ylabel('Depth relative to L4 sink (\mum)');
legend([layers + " broad", "narrow"], 'Location', 'best', 'Box', 'off');
title('Burst index vs depth');
nexttile; hold on;
for li = 1:3
    m = ok_bi & U.layer == layers(li);
    scatter(fr(m), U.bi_spont(m), 12, layer_col(li, :), 'filled', 'MarkerFaceAlpha', 0.6);
end
scatter(fr(ok_n), U.bi_spont(ok_n), 12, 'k', 'x');
yline(0, 'k:');
set(gca, 'XScale', 'log'); xlabel('Firing rate in ABCD recording (Hz)'); ylabel('BI');
r_s = corr(fr(ok_bi), U.bi_spont(ok_bi), 'type', 'Spearman');
title(sprintf('BI vs firing rate (broad rho = %.2f)', r_s));
nexttile; hold on;
ok = U.n_spont >= 100;
scatter(fr(ok & broad), U.fsb_spont(ok & broad), 12, U.bi_spont(ok & broad), 'filled');
scatter(fr(ok & narrow), U.fsb_spont(ok & narrow), 14, 'k', 'x');
rr = logspace(-1, 2, 50);
plot(rr, 1 - exp(-2 * rr * R(1).opt.isi_thr), 'b-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlabel('Firing rate (Hz)'); ylabel('FSB (spontaneous, ISI \leq 10 ms)');
colormap(gca, turbo); clim(quantile(U.bi_spont(ok_bi), [0.05 0.95])); cb = colorbar; cb.Label.String = 'BI';
legend({'broad (colour = BI)', 'narrow', 'Poisson'}, 'Location', 'northwest', 'Box', 'off');
title('FSB is rate-dependent; BI is not');
nexttile; hold on;
scatter(U.bi_spont(ok_bi), U.lv_spont(ok_bi), 12, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('BI'); ylabel('LV (spontaneous)'); yline(1, 'k:'); xline(0, 'k:');
r_s = corr(U.bi_spont(ok_bi), U.lv_spont(ok_bi), 'type', 'Spearman', 'rows', 'complete');
title(sprintf('BI vs local variation (rho = %.2f)', r_s));
nexttile; hold on;
lag = R(1).acg_t * 1000;
norm_acg = U.acg_spont ./ max(mean(U.acg_spont(:, lag >= 20 & lag < 40), 2), eps);
for b = 1:3
    plot(lag, mean(norm_acg(bicls == b, :), 1, 'omitnan'), 'Color', bi_col(b, :), 'LineWidth', 1.5);
end
plot(lag, mean(norm_acg(ok_n, :), 1, 'omitnan'), 'k--', 'LineWidth', 1.2);
yline(1, 'k:');
xlabel('Lag (ms)'); ylabel('ACG / ACG(20-40 ms)');
legend([strcat("BI ", bi_names), "narrow"], 'Box', 'off');
title('Mean spontaneous autocorrelogram');
nexttile; hold on;
edges = linspace(min(U.bi_spont(ok_bi)), max(U.bi_spont(ok_bi)), 25);
for li = 1:3
    m = ok_bi & U.layer == layers(li);
    histogram(U.bi_spont(m), edges, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', layer_col(li, :), 'LineWidth', 1.5);
end
xline(q, 'k:');
xlabel('BI'); ylabel('Fraction of units');
ct = arrayfun(@(l) sprintf('%s %.0f%%', l, 100 * mean(U.cell_type(broad & U.layer == l) == 3)), layers, 'UniformOutput', false);
mb = arrayfun(@(c) median(U.bi_spont(ok_bi & U.cell_type == c)), [1 3]);
title({'BI distribution by layer', ['cell\_type 3 (broad, slow ACG rise): ' strjoin(ct, ', ')], ...
    sprintf('median BI type1 %.2f, type3 %.2f', mb)});
legend(layers, 'Box', 'off');
title(tl, {'Spontaneous burstiness (inter-sequence grey + block pauses)', sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig4_burstiness']);

%% fig5: epoch metrics vs BI
mz = @(c, w) mean(Z(:, t >= w(1) & t < w(2), c), 2);
metrics = {
    'pre-A anticipation: ABCD z[-150,+25] - z[-600,-400] ms', mz(1, [-0.15 0.025]) - mz(1, [-0.6 -0.4])
    'A onset: ABCD z[30,80] ms', mz(1, [0.03 0.08])
    'post-onset quiescence: ABCD z[90,140] ms', mz(1, [0.09 0.14])
    sprintf('pre-B: ABCD z[%d,%d] ms', round(1000 * (dur - 0.07)), round(1000 * (dur + 0.025))), mz(1, [dur - 0.07, dur + 0.025])
    sprintf('omission: A\\_CD - E\\_CD z[%d,%d] ms', round(1000 * (dur + 0.02)), round(1000 * (dur + 0.12))), mz(2, [dur + 0.02, dur + 0.12]) - mz(3, [dur + 0.02, dur + 0.12])
    'familiar vs novel first element: A\_CD - E\_CD z[30,80] ms', mz(2, [0.03 0.08]) - mz(3, [0.03 0.08])
    };
fig = figure('Visible', 'off', 'Position', [50 50 1500 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:size(metrics, 1)
    nexttile; hold on;
    y = metrics{k, 2};
    leg = strings(1, 3);
    for li = 1:3
        m = ok_bi & U.layer == layers(li);
        scatter(U.bi_spont(m), y(m), 12, layer_col(li, :), 'filled', 'MarkerFaceAlpha', 0.6);
        [rho, p] = corr(U.bi_spont(m), y(m), 'type', 'Spearman', 'rows', 'complete');
        leg(li) = sprintf('%s rho=%.2f p=%.2g', layers(li), rho, p);
    end
    yline(0, 'k:');
    xlabel('BI'); title(metrics{k, 1}, 'FontSize', 9);
    legend(leg, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
end
title(tl, {'Per-unit response epochs vs spontaneous burst index (broad units)', sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig5_epochs_vs_BI']);

%% fig6: clusters of ABCD response shape (broad + narrow, responsive units)
k_clu = 6;
tw = t >= -0.3 & t < 4 * dur + 0.5;
Zs = smooth_t(Z, 0.010 / bin);
X = Zs(:, tw, 1);
resp = max(X(:, t(tw) >= 0), [], 2) >= 0.5 & U.layer ~= "unknown";
Xr = X(resp, :);
Xr = Xr ./ vecnorm(Xr, 2, 2);
[~, score] = pca(Xr, 'NumComponents', 10);
rng(0);
cl = kmeans(score, k_clu, 'Replicates', 50, 'MaxIter', 500);
% order clusters by latency of their mean peak
[~, pk] = max(splitapply(@(x) mean(x, 1), Xr, cl), [], 2);
[~, ord] = sort(pk);
ridx = find(resp);
fig = figure('Visible', 'off', 'Position', [50 50 1800 800]);
tl = tiledlayout(2, k_clu, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:k_clu
    mi = ridx(cl == ord(j));
    ax = nexttile(j);
    h = plot_conditions(t, Z(mi, :, :), dur, cond_col);
    title(sprintf('cluster %d (n=%d)', j, numel(mi)));
    if j == 1; ylabel('FR (z)'); legend(ax, h, strrep(conds, '_', '\_'), 'Location', 'northwest', 'Box', 'off'); end
    xlabel('Time (s)');
    nexttile(k_clu + j); hold on;
    frac_layer = arrayfun(@(l) mean(U.layer(mi) == l), layers);
    base_layer = arrayfun(@(l) mean(U.layer(resp) == l), layers);
    mb_ok = mi(ok_bi(mi));
    frac_bi = arrayfun(@(b) mean(bicls(mb_ok) == b), 1:3);
    frac_nar = mean(narrow(mi));
    b = bar(1:7, [frac_layer, frac_bi, frac_nar], 'FaceColor', 'flat');
    b.CData = [layer_col; bi_col; 0 0 0];
    plot(1:3, base_layer, 'k_', 'MarkerSize', 14, 'LineWidth', 2);
    plot(4:6, [1 1 1] / 3, 'k_', 'MarkerSize', 14, 'LineWidth', 2);
    plot(7, mean(narrow(resp)), 'k_', 'MarkerSize', 14, 'LineWidth', 2);
    set(gca, 'XTick', 1:7, 'XTickLabel', {'L23', 'L4', 'L56', 'BIlo', 'BImid', 'BIhi', 'narrow'});
    ylim([0 1]);
    title(sprintf('median BI %.2f (n=%d)', median(U.bi_spont(mb_ok)), numel(mb_ok)), 'FontSize', 9);
    if j == 1; ylabel('Fraction (bars) vs all responsive (ticks)'); end
end
title(tl, {sprintf('k-means (k=%d) on ABCD response shape (unit-normalised z-PSTH, PCA 10), %d responsive units', ...
    k_clu, nnz(resp)), sess_txt}, 'FontSize', 10);
save_png(fig, out_dir, [tag '_fig6_shape_clusters']);
end


function h = plot_conditions(t, Zm, dur, col)
hold on;
h = gobjects(1, size(Zm, 3));
for c = 1:size(Zm, 3)
    X = Zm(:, :, c);
    mu = mean(X, 1, 'omitnan');
    se = std(X, 0, 1, 'omitnan') / sqrt(max(size(X, 1), 1));
    fill([t, fliplr(t)], [mu - se, fliplr(mu + se)], col(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    h(c) = plot(t, mu, 'Color', col(c, :), 'LineWidth', 1.2);
end
finish_axes(gca, t, dur);
end


function finish_axes(ax, t, dur)
% Tight x, y from the plotted lines, then grey shading outside the sequence.
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
% Gaussian smoothing along dimension 2 (sig in bins).
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
U.session = repelem((1:numel(R))', arrayfun(@(r) numel(r.unit.cluster_id), R));
end


function cmap = diverging_map()
n = 128;
up = linspace(0, 1, n)';
cmap = [[up, up, ones(n, 1)]; [ones(n, 1), flipud(up), flipud(up)]];
end


function save_png(fig, out_dir, name)
exportgraphics(fig, fullfile(out_dir, [name '.png']), 'Resolution', 130);
close(fig);
end
