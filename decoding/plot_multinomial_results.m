function plot_multinomial_results(results_file, cfg)
% PLOT_MULTINOMIAL_RESULTS  Figures for the count-conditioned multinomial Bayesian
% decoder ('multinomial', fit_decoder) from a run_abcd_decoding results file.
%
%   plot_multinomial_results(fullfile(<analysis>, 'Decoding', 'decoding_GAVNIK_ABCD_1_results.mat'))
%
% Cross-validated posteriors are recomputed from the saved counts with the same
% folds as run_abcd_decoding (decode_with_null: rng(seed), whole sequences per fold),
% posterior = softmax of the log-likelihood (flat prior).
% Figures (saved next to the results file):
%   <base>_multinomial_bayes          per population: identity null vs observed,
%       confusion, latency sweep, mean time posterior, decoded vs true time
%   <base>_multinomial_single_trials  time posterior of example single sequences

if nargin < 2 || isempty(cfg); cfg = project_config(); end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
R = load(results_file);
R = R.results;
out_dir = fileparts(results_file);
base = sprintf('decoding_%s', R.stimulus);
tag = sprintf('%s %s %s', R.subject, R.date, strrep(R.stimulus, '_', ' '));
trials = R.trials;
tc = R.time_centres;
n_el = numel(trials.element_names);
n_tb = numel(tc);
pops = R.pop;
np = numel(pops);

% cross-validated posteriors
for p = 1:np
    n_seq = size(pops(p).counts_time, 1);
    Xe = reshape(permute(pops(p).counts_element, [3 1 2]), n_el * n_seq, []);
    ye = repmat((1:n_el)', n_seq, 1);
    % 5 classes A B C D grey, sequence-major as in run_abcd_decoding
    Xg = reshape(permute(cat(3, pops(p).counts_element, pops(p).counts_grey), [3 1 2]), (n_el + 1) * n_seq, []);
    yg = repmat((1:n_el + 1)', n_seq, 1);
    post(p).grey = cv_posterior(Xg, yg, repelem((1:n_seq)', n_el + 1), R.opts); %#ok<AGROW>
    Xt = reshape(permute(pops(p).counts_time, [3 1 2]), n_tb * n_seq, []);
    yt = repmat((1:n_tb)', n_seq, 1);
    post(p).el = cv_posterior(Xe, ye, repelem((1:n_seq)', n_el), R.opts); %#ok<AGROW>
    post(p).time = cv_posterior(Xt, yt, repelem((1:n_seq)', n_tb), R.opts); %#ok<AGROW>
    post(p).yt = yt; %#ok<AGROW>
    [~, pr] = max(post(p).el, [], 2);
    m = strcmp({pops(p).identity.method}, 'multinomial');
    fprintf('%s: recomputed identity accuracy %.3f (saved %.3f)\n', pops(p).name, mean(pr == ye), pops(p).identity(m).acc);
end

figs = [plot_summary(pops, post, trials, tc, R.opts, tag), plot_single_trials(pops, post, trials, tc, tag)];
names = {[base '_multinomial_bayes'], [base '_multinomial_single_trials']};
save_vector_figures(figs, names, out_dir, cfg);
end


function P = cv_posterior(X, y, group, opts)
% same folds as decode_with_null
rng(opts.seed);
[~, ~, g] = unique(group);
gf = mod(randperm(max(g)) - 1, opts.n_folds) + 1;
fold = gf(g)';
K = max(y);
P = zeros(numel(y), K);
for f = 1:opts.n_folds
    test = fold == f;
    m = fit_decoder(X(~test, :), y(~test), K, 'multinomial');
    [~, ll] = apply_decoder(m, X(test, :));
    ll = ll - max(ll, [], 2);
    P(test, :) = exp(ll) ./ sum(exp(ll), 2);
end
end


%% ------------------------------------------------------------------------
function cmap = seq_cmap()
cmap = interp1([0 0.6 1], [1 1 1; 0.165 0.471 0.839; 0.04 0.15 0.35], linspace(0, 1, 256));
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function mark_sequence(ax, trials, vertical_too)
s = 1000 * unique([trials.element_start trials.element_end]);
xline(ax, s, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
if vertical_too; yline(ax, s, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); end
end


function f = plot_summary(pops, post, trials, tc, opts, tag)
np = numel(pops);
n_el = numel(trials.element_names);
f = figure('Color', 'w', 'Position', [40 40 1700 330 * np]);
tl = tiledlayout(f, np, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - multinomial Bayesian decoder (conditioned on total spike count: gain-invariant), ' ...
    '%d-fold CV over sequences'], tag, opts.n_folds), 'FontSize', 11);
blue = [0.165 0.471 0.839]; grey = [0.55 0.55 0.55];
for p = 1:np
    P = pops(p);
    mi = strcmp({P.identity.method}, 'multinomial');
    bi = strcmp({P.identity.method}, 'bayes');
    pname = sprintf('%s (n = %d)', strrep(P.name, '_', ' '), numel(P.cluster_id));

    % 1. identity: observed vs null
    ax = nexttile(tl); hold(ax, 'on');
    edges = 0.15:0.005:max(0.45, P.identity(mi).acc + 0.03);
    histogram(ax, 100 * P.identity(mi).null_acc, 100 * edges, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', ...
        'DisplayName', 'null (labels shuffled within sequence)');
    yl = ylim(ax);
    xline(ax, 100 * P.identity(mi).acc, '-', 'Color', blue, 'LineWidth', 2.5, 'DisplayName', ...
        sprintf('stimulus %.1f%% (p = %.3f)', 100 * P.identity(mi).acc, P.identity(mi).p));
    xline(ax, 100 * P.identity_pre(mi).acc, '-', 'Color', grey, 'LineWidth', 2.5, 'DisplayName', ...
        sprintf('grey control %.1f%% (p = %.3f)', 100 * P.identity_pre(mi).acc, P.identity_pre(mi).p));
    xline(ax, 100 * P.identity(bi).acc, '--', 'Color', blue, 'LineWidth', 1.2, 'DisplayName', ...
        sprintf('Poisson Bayes, stimulus %.1f%%', 100 * P.identity(bi).acc));
    xline(ax, 25, ':', 'Color', 'k', 'HandleVisibility', 'off');
    ylim(ax, yl .* [1 1.9]);
    if P.identity(mi).acc > 0.6; xlim(ax, [15 100]); end
    xlabel(ax, 'A/B/C/D accuracy (%), chance 25'); ylabel(ax, 'null count');
    legend(ax, 'Location', 'north', 'Box', 'off', 'FontSize', 6.5);
    title(ax, [pname ': element identity'], 'FontWeight', 'normal');
    style(ax);

    % 2. A/B/C/D/grey: mean posterior per true class
    ax = nexttile(tl);
    nc = n_el + 1;
    yg = repmat((1:nc)', size(post(p).grey, 1) / nc, 1);
    M = zeros(nc);
    for k = 1:nc
        M(k, :) = mean(post(p).grey(yg == k, :), 1);
    end
    imagesc(ax, 100 * M, [0 max(50, 100 * max(M(:)))]); colormap(ax, seq_cmap()); axis(ax, 'square');
    for i = 1:nc
        for j = 1:nc
            text(ax, j, i, sprintf('%.0f', 100 * M(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', [1 1 1] * (M(i, j) > 0.55));
        end
    end
    cls = [trials.element_names, {'grey'}];
    set(ax, 'XTick', 1:nc, 'XTickLabel', cls, 'YTick', 1:nc, 'YTickLabel', cls);
    xlabel(ax, 'decoded'); ylabel(ax, 'true');
    cb = colorbar(ax); cb.Label.String = 'mean posterior (%)';
    gi = strcmp({P.identity_grey.method}, 'multinomial');
    title(ax, sprintf('A/B/C/D/grey mean posterior: %.0f%% correct (chance 20%%, p = %.3f)', ...
        100 * P.identity_grey(gi).acc, P.identity_grey(gi).p), 'FontWeight', 'normal');
    style(ax);

    % 3. latency sweep, multinomial vs Poisson Bayes
    ax = nexttile(tl); hold(ax, 'on');
    S = P.latency_sweep;
    lat = 1000 * S.latency;
    km = find(strcmp(S.methods, 'multinomial')); kb = find(strcmp(S.methods, 'bayes'));
    plot(ax, lat, 100 * S.null95(km, :), '--', 'Color', blue, 'HandleVisibility', 'off');
    plot(ax, lat, 100 * S.acc(kb, :), '-', 'Color', grey, 'LineWidth', 1.5, 'DisplayName', 'Poisson Bayes');
    plot(ax, lat, 100 * S.acc(km, :), '-o', 'Color', blue, 'MarkerFaceColor', blue, 'MarkerSize', 4, ...
        'LineWidth', 2, 'DisplayName', 'multinomial');
    xline(ax, 0, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    yline(ax, 25, ':', 'Color', 'k', 'HandleVisibility', 'off');
    ylim(ax, [0 100]);
    xlabel(ax, 'window start re element onset (ms)'); ylabel(ax, 'accuracy (%)');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, 'latency sweep (150 ms window; dashed = null 95%)', 'FontWeight', 'normal');
    style(ax);

    % 4. time: mean posterior, true x decoded
    ax = nexttile(tl);
    n_tb = numel(tc);
    T = zeros(n_tb);
    for k = 1:n_tb
        T(k, :) = mean(post(p).time(post(p).yt == k, :), 1);
    end
    imagesc(ax, 1000 * tc, 1000 * tc, 100 * T, [0 max(15, 100 * prctile(T(:), 99))]);
    colormap(ax, seq_cmap()); axis(ax, 'square'); set(ax, 'YDir', 'normal'); hold(ax, 'on');
    plot(ax, 1000 * tc([1 end]), 1000 * tc([1 end]), ':', 'Color', [0.4 0.4 0.4]);
    mark_sequence(ax, trials, true);
    cb = colorbar(ax); cb.Label.String = 'mean posterior (%)';
    mt = strcmp({P.time.method}, 'multinomial');
    xlabel(ax, 'decoded time (ms)'); ylabel(ax, 'true time from A onset (ms)');
    title(ax, sprintf('time in sequence: %.1f%% exact (chance %.0f%%)', 100 * P.time(mt).acc, 100 / n_tb), ...
        'FontWeight', 'normal');
    style(ax);

    % 5. decoded (posterior mean) vs true time
    ax = nexttile(tl); hold(ax, 'on');
    pm = post(p).time * tc(:);                         % posterior-mean decoded time per sample
    q = zeros(n_tb, 3);
    for k = 1:n_tb
        q(k, :) = prctile(pm(post(p).yt == k), [25 50 75]);
    end
    fill(ax, 1000 * [tc fliplr(tc)], 1000 * [q(:, 1)' fliplr(q(:, 3)')], blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'DisplayName', 'IQR over sequences');
    plot(ax, 1000 * tc, 1000 * q(:, 2), '-o', 'Color', blue, 'MarkerFaceColor', blue, 'MarkerSize', 4, ...
        'LineWidth', 2, 'DisplayName', 'median');
    plot(ax, 1000 * tc([1 end]), 1000 * tc([1 end]), 'k:', 'DisplayName', 'perfect');
    yline(ax, 1000 * mean(tc), '--', 'Color', grey, 'DisplayName', 'no information');
    mark_sequence(ax, trials, false);
    axis(ax, 'square'); xlim(ax, 1000 * tc([1 end])); ylim(ax, 1000 * tc([1 end]));
    xlabel(ax, 'true time (ms)'); ylabel(ax, 'decoded time, posterior mean (ms)');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('decoded vs true, MAE %.0f ms', 1000 * P.time(mt).mae), 'FontWeight', 'normal');
    style(ax);
end
end


function f = plot_single_trials(pops, post, trials, tc, tag)
% Time posterior (decoded time x true time) of single sequences: what one
% presentation looks like to the decoder.
np = numel(pops);
n_ex = 8;
n_tb = numel(tc);
n_seq = size(post(1).time, 1) / n_tb;
ex = round(linspace(1, n_seq, n_ex));
f = figure('Color', 'w', 'Position', [40 40 1700 260 * np]);
tl = tiledlayout(f, np, n_ex, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - multinomial decoder, single sequences: posterior over decoded time (y) at each true time (x); ' ...
    'diagonal = perfect'], tag), 'FontSize', 11);
for p = 1:np
    for i = 1:n_ex
        ax = nexttile(tl); hold(ax, 'on');
        rows = (ex(i) - 1) * n_tb + (1:n_tb);
        imagesc(ax, 1000 * tc, 1000 * tc, post(p).time(rows, :)', [0 1]);
        colormap(ax, seq_cmap()); set(ax, 'YDir', 'normal'); axis(ax, 'square', 'tight');
        plot(ax, 1000 * tc([1 end]), 1000 * tc([1 end]), ':', 'Color', [0.5 0.5 0.5]);
        mark_sequence(ax, trials, true);
        title(ax, sprintf('sequence %d', ex(i)), 'FontWeight', 'normal', 'FontSize', 8);
        if i == 1
            ylabel(ax, {strrep(pops(p).name, '_', ' '), 'decoded time (ms)'});
        end
        if p == np; xlabel(ax, 'true time (ms)'); end
        style(ax);
    end
end
cb = colorbar(ax); cb.Label.String = 'posterior';
end
