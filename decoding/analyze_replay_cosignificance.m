function R = analyze_replay_cosignificance(subject, date, rec_name, stimulus_name, cfg, opts)
% ANALYZE_REPLAY_COSIGNIFICANCE  Lseq replay: how many events are significant in CA1 and in V1,
% whether CA1 and V1 are significant in the SAME events more often than expected, and example
% rasters of significant events.
%
%   R = analyze_replay_cosignificance('M00069', '20250929', 'Sleep_Box_2', 'GAVNIK_ABCD_1', project_config())
%   R = analyze_replay_cosignificance('M00069', '20250929', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_1', ...)  % awake
%
% Reads run_sleep_replay results (sleep recordings) or run_abcd_awake_replay results (when
% REC_NAME is the stimulus recording). CA1 (CA1_all) and V1 (V1_good) are scored on the same
% events, each with its own awake order matrix.
%
% Significant = calibrated (cell-ID shuffle AND beyond the matched control windows), either
% direction, alpha 0.025 per tail as in the replay scripts; 'raw' = cell-ID shuffle only.
% Co-significance, per event group (sleep: ripples / merged clusters / MUA bursts without
% ripple, in opts.states; awake: all ripple-MUA candidates), events scored in both regions:
%   observed % both significant vs expected under independence (pCA1 x pV1), odds ratio and
%   one-sided Fisher test; and a matched null: V1 labels permuted among events of similar size
%   (strata = terciles of duration x CA1 spike count x V1 spike count), because longer events
%   with more spikes can be significant in both regions for that reason alone.
%   Also the rank correlation of CA1 z and V1 z, plain and partial (controlling for log
%   duration and log spike counts), which uses every event rather than the few significant ones.
% Examples: up to 12 events (co-significant first, then CA1-only, then V1-only), CA1 and V1
%   rasters with neurons ordered by the awake order index and spikes coloured by it; ripple
%   windows shaded.
%
% Outputs: <rec>\replay\replay_<rec>_lseq_cosig.mat / .fig / .pdf and replay_<rec>_lseq_examples

if nargin < 5 || isempty(cfg); cfg = project_config(); end
if nargin < 6; opts = struct(); end
def.states = {'SWS'};
def.n_perm = 5000;
def.n_examples = 12;
def.pad = 0.05;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
awake = strcmp(rec_name, stimulus_name);
paths = get_session_paths(subject, date, rec_name, cfg);
out_dir = fullfile(paths.analysis_dir, 'replay');
check_output_path(out_dir, cfg);

%% events and scores
if awake
    X = load(fullfile(out_dir, sprintf('replay_%s_results.mat', rec_name))).results;
    ev = X.scores(1).event;
    E = table();
    E.onset = X.candidates.onset(ev); E.offset = X.candidates.offset(ev);
    E.duration = E.offset - E.onset;
    E.kind = repmat("candidate", numel(ev), 1);
    E.state = repmat("awake", numel(ev), 1);
    rip_iv = [X.ripples.onset(:) X.ripples.offset(:)];
    groups = {"all candidates", true(numel(ev), 1)};
    alpha = X.opts.alpha;
else
    X = load(fullfile(out_dir, sprintf('replay_%s_template_%s_results.mat', rec_name, stimulus_name))).results;
    ev = X.scores(1).event;
    E = X.events(ev, :);
    rip_iv = X.clusters.window;
    in_state = ismember(E.state, string(opts.states));
    groups = {"ripples", in_state & E.kind == "ripple"; "merged clusters", in_state & E.kind == "cluster"; ...
        "MUA bursts, no ripple", in_state & E.kind == "burst"};
    alpha = X.opts.alpha;
end
iC = find(strcmp({X.scores.name}, 'CA1_all')); iV = find(strcmp({X.scores.name}, 'V1_good'));
assert(isequal(X.scores(iC).event, X.scores(iV).event), 'CA1 and V1 were scored on different events');
sC = X.scores(iC); sV = X.scores(iV);
sigC = sC.forward_cal | sC.reverse_cal; sigV = sV.forward_cal | sV.reverse_cal;
rawC = sC.forward | sC.reverse;         rawV = sV.forward | sV.reverse;
dirC = sign(sC.forward_cal - sC.reverse_cal); dirV = sign(sV.forward_cal - sV.reverse_cal);
both = sC.scored & sV.scored;
tag = sprintf('%s %s %s (template %s)', subject, date, strrep(rec_name, '_', ' '), strrep(stimulus_name, '_', ' '));

%% proportions and co-significance per group
rng(opts.seed);
G = table();
for g = 1:size(groups, 1)
    m = groups{g, 2};
    r = struct();
    r.group = string(groups{g, 1});
    r.n_CA1 = sum(m & sC.scored); r.n_V1 = sum(m & sV.scored);
    r.CA1_fwd = mean(sC.forward_cal(m & sC.scored)); r.CA1_rev = mean(sC.reverse_cal(m & sC.scored));
    r.V1_fwd = mean(sV.forward_cal(m & sV.scored));  r.V1_rev = mean(sV.reverse_cal(m & sV.scored));
    r.CA1_raw = mean(rawC(m & sC.scored)); r.V1_raw = mean(rawV(m & sV.scored));
    b = m & both;
    r.n_both_scored = sum(b);
    [r.pCA1, r.pV1, r.p_joint, r.p_expected, r.odds_ratio, r.p_fisher, r.p_joint_null, r.p_joint_null95, r.p_matched, ...
        r.n_cosig, r.n_same_dir] = cosig(sigC(b), sigV(b), dirC(b), dirV(b), E.duration(b), sC.n_spikes(b), sV.n_spikes(b), opts.n_perm);
    [~, ~, r.raw_p_joint, r.raw_p_expected, r.raw_odds_ratio, r.raw_p_fisher, ~, ~, r.raw_p_matched] = ...
        cosig(rawC(b), rawV(b), [], [], E.duration(b), sC.n_spikes(b), sV.n_spikes(b), opts.n_perm);
    [r.rho_z, r.p_rho_z, r.partial_rho_z, r.p_partial_rho_z] = deal(NaN);
    if sum(b) > 10
        [r.rho_z, r.p_rho_z] = corr(sC.z(b), sV.z(b), 'Type', 'Spearman');
        cov_ = [log(E.duration(b)), log(sC.n_spikes(b) + 1), log(sV.n_spikes(b) + 1)];
        [r.partial_rho_z, r.p_partial_rho_z] = partialcorr(sC.z(b), sV.z(b), cov_, 'Type', 'Spearman');
    end
    G = [G; struct2table(r)]; %#ok<AGROW>
end
disp(G(:, {'group', 'n_both_scored', 'pCA1', 'pV1', 'p_joint', 'p_expected', 'p_joint_null', 'n_cosig', 'p_fisher', 'p_matched', ...
    'raw_p_joint', 'raw_p_expected', 'raw_p_matched', 'partial_rho_z', 'p_partial_rho_z'}));

R = struct('subject', subject, 'date', date, 'recording', rec_name, 'stimulus', stimulus_name, 'awake', awake, ...
    'opts', opts, 'alpha', alpha, 'groups', G, 'events', E, 'zCA1', sC.z, 'zV1', sV.z, 'sigCA1', sigC, 'sigV1', sigV, ...
    'dirCA1', dirC, 'dirV1', dirV, 'scored_both', both, 'nspkCA1', sC.n_spikes, 'nspkV1', sV.n_spikes);

%% figures
f1 = plot_summary(R, groups, tag);
[f2, R.example_events] = plot_examples(subject, date, rec_name, cfg, X, sC, sV, E, rip_iv, sigC, sigV, dirC, dirV, ...
    both & any([groups{:, 2}], 2), opts, tag, awake);
if ~opts.save; return; end
base = sprintf('replay_%s_lseq', rec_name);
safe_save(fullfile(out_dir, [base '_cosig.mat']), struct('results', R), cfg);
save_vector_figures([f1 f2], {[base '_cosig'], [base '_examples']}, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function [pC, pV, pJ, pE, orr, pF, pNull, pNull95, pMatched, nJ, nSame] = cosig(a, b, da, db, dur, nC, nV, n_perm)
[pC, pV, pJ, pE, orr, pF, pNull, pNull95, pMatched, nJ, nSame] = deal(NaN);
if numel(a) < 10; return; end
pC = mean(a); pV = mean(b); pJ = mean(a & b); pE = pC * pV; nJ = sum(a & b);
if ~isempty(da); nSame = sum(a & b & da == db); end
t = [sum(a & b) sum(a & ~b); sum(~a & b) sum(~a & ~b)];
[~, pF, st] = fishertest(t, 'Tail', 'right'); orr = st.OddsRatio;
% matched null: permute V1 labels within strata of similar event size
q = @(x) discretize(tiedrank(x) / numel(x), [0 1/3 2/3 1]);
[~, ~, strat] = unique([q(dur) q(nC) q(nV)], 'rows');
null = zeros(n_perm, 1);
for k = 1:n_perm
    bp = b;
    for s = 1:max(strat)
        i = find(strat == s);
        bp(i) = b(i(randperm(numel(i))));
    end
    null(k) = mean(a & bp);
end
pNull = mean(null); pNull95 = prctile(null, 95);
pMatched = (1 + sum(null >= pJ)) / (1 + n_perm);
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function [lo, hi] = wilson(p, n)
z = 1.96;
c = (p + z^2 ./ (2 * n)) ./ (1 + z^2 ./ n);
h = z * sqrt(p .* (1 - p) ./ n + z^2 ./ (4 * n.^2)) ./ (1 + z^2 ./ n);
lo = c - h; hi = c + h;
end

function f = plot_summary(R, groups, tag)
G = R.groups; ng = height(G);
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.55 0.55 0.55]; green = [0.106 0.686 0.478];
f = figure('Color', 'w', 'Position', [40 40 1700 820]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - Lseq replay: proportion significant and CA1-V1 co-significance (calibrated, alpha %.3g per tail)', ...
    tag, R.alpha), 'FontSize', 10.5);

% 1. proportions significant
ax = nexttile(tl, [1 2]); hold(ax, 'on');
V = 100 * [G.CA1_fwd G.CA1_rev G.V1_fwd G.V1_rev];
N = [G.n_CA1 G.n_CA1 G.n_V1 G.n_V1];
bb = bar(ax, 1:ng, V, 'BarWidth', 0.85, 'EdgeColor', 'none');
cols = [blue; blue * 0.4 + 0.6; orange; orange * 0.4 + 0.6];
for k = 1:4
    bb(k).FaceColor = cols(k, :);
    [lo, hi] = wilson(V(:, k) / 100, N(:, k));
    errorbar(ax, bb(k).XEndPoints, V(:, k), V(:, k) - 100 * lo, 100 * hi - V(:, k), 'k', 'LineStyle', 'none', ...
        'CapSize', 3, 'HandleVisibility', 'off');
    text(ax, bb(k).XEndPoints, 100 * hi(:)' + 0.4, compose('%d', round(V(:, k)' .* N(:, k)' / 100)), ...
        'HorizontalAlignment', 'center', 'FontSize', 6.5, 'Color', grey);
end
yline(ax, 100 * R.alpha, ':', sprintf('%.1f%% per tail', 100 * R.alpha), 'FontSize', 7, 'HandleVisibility', 'off');
set(ax, 'XTick', 1:ng, 'XTickLabel', compose('%s (n = %d / %d)', G.group, G.n_CA1, G.n_V1));
ylabel(ax, '% of scored events significant');
legend(ax, {'CA1 forward', 'CA1 reverse', 'V1 forward', 'V1 reverse'}, 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 7);
title(ax, 'calibrated significance (numbers = events); n = CA1 / V1 events scored', 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

% 2. co-significance
ax = nexttile(tl); hold(ax, 'on');
O = 100 * [G.p_joint G.p_expected G.p_joint_null];
bb = bar(ax, 1:ng, O, 'BarWidth', 0.85, 'EdgeColor', 'none');
bb(1).FaceColor = green; bb(2).FaceColor = grey; bb(3).FaceColor = [0.8 0.8 0.8];
errorbar(ax, bb(3).XEndPoints, O(:, 3), zeros(ng, 1), 100 * G.p_joint_null95 - O(:, 3), 'k', 'LineStyle', 'none', 'CapSize', 3);
for g = 1:ng
    text(ax, bb(1).XEndPoints(g), max(O(g, :)) * 1.08 + 0.05, sprintf('%d ev; OR %.1f\nFisher p %.2g\nmatched p %.2g', ...
        G.n_cosig(g), G.odds_ratio(g), G.p_fisher(g), G.p_matched(g)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 6.5);
end
ylim(ax, [0 max([O(:); 0.5]) * 1.6]);
set(ax, 'XTick', 1:ng, 'XTickLabel', G.group); xtickangle(ax, 15);
ylabel(ax, '% of events significant in BOTH');
legend(ax, {'observed', 'expected (pCA1 x pV1)', 'matched-size null (mean, 95%)'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
title(ax, 'CA1 and V1 significant in the same event', 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

% 3. z scatter per group
for g = 1:min(ng, 3)
    ax = nexttile(tl); hold(ax, 'on');
    m = groups{g, 2} & R.scored_both;
    zc = R.zCA1(m); zv = R.zV1(m); sc = R.sigCA1(m); sv = R.sigV1(m);
    scatter(ax, zc(~sc & ~sv), zv(~sc & ~sv), 6, grey, 'filled', 'MarkerFaceAlpha', 0.4);
    scatter(ax, zc(sc & ~sv), zv(sc & ~sv), 14, blue, 'filled');
    scatter(ax, zc(~sc & sv), zv(~sc & sv), 14, orange, 'filled');
    scatter(ax, zc(sc & sv), zv(sc & sv), 30, green, 'filled', 'MarkerEdgeColor', 'k');
    xline(ax, 0, 'k:'); yline(ax, 0, 'k:');
    xlabel(ax, 'CA1 Lseq z'); ylabel(ax, 'V1 Lseq z');
    title(ax, sprintf('%s: rho %.2f (p %.2g); partial rho %.2f (p %.2g)', G.group(g), G.rho_z(g), G.p_rho_z(g), ...
        G.partial_rho_z(g), G.p_partial_rho_z(g)), 'FontWeight', 'normal', 'FontSize', 8);
    if g == 1
        legend(ax, {'neither', 'CA1 only', 'V1 only', 'both'}, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
    end
    style(ax);
end
end

function [f, picked] = plot_examples(subject, date, rec_name, cfg, X, sC, sV, E, rip_iv, sigC, sigV, dirC, dirV, pool, opts, tag, awake)
% pick events: co-significant, CA1 only, V1 only (largest |z| first)
zC = abs(sC.z); zV = abs(sV.z);
cand = {find(pool & sigC & sigV), min(zC, zV); find(pool & sigC & ~sigV), zC; find(pool & ~sigC & sigV), zV};
picked = []; lab = strings(0);
names = ["CA1 + V1", "CA1 only", "V1 only"];
per = ceil(opts.n_examples / 3);
for c = 1:3
    i = cand{c, 1}; [~, o] = sort(cand{c, 2}(i), 'descend'); i = i(o);
    take = i(1:min(per, numel(i)));
    picked = [picked; take(:)]; lab = [lab; repmat(names(c), numel(take), 1)]; %#ok<AGROW>
end
spare = setdiff([cand{1, 1}; cand{2, 1}; cand{3, 1}], picked, 'stable');
while numel(picked) < opts.n_examples && ~isempty(spare)
    e = spare(1); spare(1) = [];
    picked(end + 1, 1) = e; %#ok<AGROW>
    if sigC(e) && sigV(e); lab(end + 1, 1) = names(1); elseif sigC(e); lab(end + 1, 1) = names(2); else; lab(end + 1, 1) = names(3); end
end
f = figure('Color', 'w', 'Position', [20 20 1800 1000]);
ncol = 4; nrow = max(1, ceil(numel(picked) / ncol));
tl = tiledlayout(f, 2 * nrow, ncol, 'TileSpacing', 'tight', 'Padding', 'compact');
title(tl, sprintf('%s - example Lseq events (rows: active neurons sorted by awake order; colour = awake order rank, blue early -> red late; grey = ripple windows; dotted = scored window)', tag), ...
    'FontSize', 10.5);
if isempty(picked)
    ax = nexttile(tl); axis(ax, 'off'); text(ax, 0.5, 0.5, 'no significant events', 'HorizontalAlignment', 'center');
    return
end
paths = get_session_paths(subject, date, rec_name, cfg);
units = load_session_clusters(paths, get_probe_layout(paths, cfg), cfg);
cmap = turbo(256);
dirs = ["rev", "", "fwd"];
for k = 1:numel(picked)
    e = picked(k);
    t0 = E.onset(e) - opts.pad; t1 = E.offset(e) + opts.pad;
    for r = 1:2
        if r == 1; s = sC; rname = 'CA1'; else; s = sV; rname = 'V1'; end
        w = collect_spike_windows(units.spike_times, units.spike_id, s.cluster_id, t0, t1); w = w{1};
        % rows: neurons that fire in the window, ordered by the awake order index; colour = that
        % neuron's rank in the awake order among all neurons (early = blue, late = red)
        n_all = numel(s.order_index);
        rank_all = tiedrank(s.order_index(:));
        act = unique(w(:, 2));
        [~, o] = sort(s.order_index(act)); row = zeros(n_all, 1); row(act(o)) = 1:numel(act);
        n_row = max(numel(act), 1);
        ax = nexttile(tl, (2 * floor((k - 1) / ncol) + r - 1) * ncol + mod(k - 1, ncol) + 1); hold(ax, 'on');
        ri = find(rip_iv(:, 2) > t0 & rip_iv(:, 1) < t1);
        for j = ri'
            x = 1000 * ([max(rip_iv(j, 1), t0) min(rip_iv(j, 2), t1)] - t0);
            patch(ax, x([1 2 2 1]), [0 0 n_row + 1 n_row + 1], [0.92 0.92 0.92], 'EdgeColor', 'none');
        end
        xline(ax, 1000 * opts.pad, 'k:'); xline(ax, 1000 * (E.offset(e) - t0), 'k:');
        if ~isempty(w)
            ci = max(1, min(256, round(1 + 255 * (rank_all(w(:, 2)) - 1) / max(n_all - 1, 1))));
            scatter(ax, 1000 * w(:, 1), row(w(:, 2)), 12, cmap(ci, :), '|', 'LineWidth', 1.2);
        end
        xlim(ax, [0 1000 * (t1 - t0)]); ylim(ax, [0 n_row + 1]);
        set(ax, 'YTick', []); style(ax);
        ylabel(ax, rname, 'FontSize', 7);
        if r == 1
            title(ax, sprintf('%s | %s %s %.0f ms\nCA1 z %.1f%s%s | V1 z %.1f%s%s', lab(k), E.kind(e), E.state(e), ...
                1000 * E.duration(e), sC.z(e), dirs(dirC(e) + 2), repmat('*', 1, sigC(e)), sV.z(e), dirs(dirV(e) + 2), ...
                repmat('*', 1, sigV(e))), 'FontWeight', 'normal', 'FontSize', 7.5);
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'ms from window start', 'FontSize', 7);
        end
    end
end
end
