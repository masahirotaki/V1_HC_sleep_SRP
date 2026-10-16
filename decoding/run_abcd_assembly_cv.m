function results = run_abcd_assembly_cv(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_ASSEMBLY_CV  Can cell-assembly activation (PCA -> Marchenko-Pastur -> ICA;
% Lopes-dos-Santos et al. 2013) learned on some awake sequences tell A / B / C / D / grey
% apart on other awake sequences? This is the awake test of the sleep plan: detect assembly
% reactivation per 20 ms bin, then say which element it was.
%
%   results = run_abcd_assembly_cv('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Data: spike counts in opts.bin (20 ms) bins from -0.6 to +0.8 s around each A onset.
% Labelled bins: A-D = bins whose centre lies in the element's on-time; grey = -0.6 to
% -0.45 s before A (one 150 ms window per sequence, as the 5-class decoding). Windows: the
% five 150 ms windows per sequence (mean activation over their bins).
% 10-fold CV over sequences. On the training sequences: z-score each neuron (all bins of the
% window), find assemblies, label them by their mean activation in A-D / grey. Held-out
% sequences are scored with:
%   ica_trainz     ICA templates; test bins z-scored with the training mean / SD
%   ica_ownz       test bins z-scored with their own mean / SD (as sleep: own-epoch z-score)
%   ica_ownz_cm    as ica_ownz, plus the common mode (population mean z per bin) removed
%                  from training and test bins (against shared gain changes, e.g. ripples)
%   pca_ownz       PCA templates (significant PCs) instead of ICA, own z-score
%   multinomial    no assemblies: count-conditioned multinomial Bayes on the raw counts
% Classification: shrinkage LDA on the assembly activation vector (fit_decoder 'lda'),
% per single bin and per 150 ms window; chance 20%. Null (windows): training labels
% permuted within each sequence, opts.n_perm times.
% Outputs (decoding_* names) in <analysis>/<subject>/<date>/<stimulus>/Decoding.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.populations = struct('name', {'CA1_all', 'V1_good'}, 'region', {"CA1", "V1"}, 'good_only', {false, true});
def.bin = 0.02;
def.range = [];               % [] = [-0.6, sequence + 0.2] s around A onset (-0.6 to 0.8 s for 150 ms)
def.grey = [-0.6 -0.45];
def.n_folds = 10;
def.n_perm = 100;
def.min_spikes = 20;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths  = get_session_paths(subject, date, stimulus_name, cfg);
units  = load_session_clusters(paths, get_probe_layout(paths, cfg), cfg);
trials = get_abcd_trials(paths);
out_dir = fullfile(paths.analysis_dir, 'Decoding');
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));
n_seq = numel(trials.seq_onset);
n_el = numel(trials.element_names);
cls = [trials.element_names, {'grey'}];
if isempty(opts.range); opts.range = [-0.6, round(trials.seq_duration, 2) + 0.2]; end
nc = n_el + 1;

edges = opts.range(1):opts.bin:opts.range(2) + 1e-9;
tc = edges(1:end - 1) + opts.bin / 2;
nb = numel(tc);
lab = zeros(1, nb);
for e = 1:n_el
    lab(tc >= trials.element_start(e) & tc < trials.element_start(e) + trials.element_duration) = e;
end
lab(tc >= opts.grey(1) & tc < opts.grey(2)) = nc;
rng(opts.seed);
fold = mod(randperm(n_seq) - 1, opts.n_folds) + 1;
variants = {'ica_trainz', 'ica_ownz', 'ica_ownz_cm', 'pca_ownz', 'multinomial'};
vlab = {'ICA, train z', 'ICA, own z', 'ICA, own z, no common mode', 'PCA, own z', 'multinomial (counts)'};

for p = 1:numel(opts.populations)
    P = opts.populations(p);
    ids = units.cluster_id(units.region == P.region & (units.is_good | ~P.good_only));
    X = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset, edges);   % seq x neuron x bin
    keep = squeeze(sum(X, [1 3])) >= opts.min_spikes;
    X = X(:, keep, :); ids = ids(keep);
    n = numel(ids);
    fprintf('\n=== %s: %d neurons, %d sequences, %d bins of %.0f ms\n', P.name, n, n_seq, nb, 1000 * opts.bin);
    Xb = reshape(permute(X, [3 1 2]), nb * n_seq, n);          % bins (sequence-major) x neuron
    seq_of_bin = repelem((1:n_seq)', nb);
    lab_of_bin = repmat(lab', n_seq, 1);

    res = struct();
    for v = 1:numel(variants)
        res.(variants{v}) = struct('bin_pred', zeros(nb * n_seq, 1), 'win_pred', zeros(n_seq * nc, 1));
    end
    win_true = repmat((1:nc)', n_seq, 1);
    win_seq = repelem((1:n_seq)', nc);
    n_asm = zeros(opts.n_folds, 1);
    perm_acc = zeros(opts.n_perm, 2);       % ica_ownz, multinomial
    perm_pred = cell(opts.n_perm, 2);
    for f = 1:opts.n_folds
        tr = fold(seq_of_bin) ~= f;
        te = ~tr;
        mu = mean(Xb(tr, :), 1); sd = std(Xb(tr, :), 0, 1); sd(sd == 0) = 1;
        Ztr = (Xb(tr, :) - mu) ./ sd;
        mu2 = mean(Xb(te, :), 1); sd2 = std(Xb(te, :), 0, 1); sd2(sd2 == 0) = 1;
        Zte_train = (Xb(te, :) - mu) ./ sd;
        Zte_own = (Xb(te, :) - mu2) ./ sd2;
        A_ica = assembly_templates(Ztr, 'ica', opts.seed);
        A_cm = assembly_templates(Ztr - mean(Ztr, 2), 'ica', opts.seed);
        A_pca = assembly_templates(Ztr, 'pca');
        n_asm(f) = size(A_ica.W, 2);
        F.ica_trainz = {assembly_activation(Ztr, A_ica.W), assembly_activation(Zte_train, A_ica.W)};
        F.ica_ownz = {assembly_activation(Ztr, A_ica.W), assembly_activation(Zte_own, A_ica.W)};
        F.ica_ownz_cm = {assembly_activation(Ztr - mean(Ztr, 2), A_cm.W), assembly_activation(Zte_own - mean(Zte_own, 2), A_cm.W)};
        F.pca_ownz = {assembly_activation(Ztr, A_pca.W), assembly_activation(Zte_own, A_pca.W)};
        F.multinomial = {Xb(tr, :), Xb(te, :)};
        ltr = lab_of_bin(tr); lte = lab_of_bin(te);
        str = seq_of_bin(tr); ste = seq_of_bin(te);
        for v = 1:numel(variants)
            vn = variants{v};
            Ftr = F.(vn){1}; Fte = F.(vn){2};
            if isempty(Ftr); continue; end
            method = 'lda'; if strcmp(vn, 'multinomial'); method = 'multinomial'; end
            % single bins
            m = fit_decoder(Ftr(ltr > 0, :), ltr(ltr > 0), nc, method);
            pr = zeros(sum(te), 1);
            pr(lte > 0) = apply_decoder(m, Fte(lte > 0, :));
            idx = find(te);
            res.(vn).bin_pred(idx) = pr;
            % 150 ms windows: mean activation (or summed counts) over the window's bins
            [Wtr, ytr] = windows(Ftr, ltr, str, nc, strcmp(vn, 'multinomial'));
            [Wte, yte, sq] = windows(Fte, lte, ste, nc, strcmp(vn, 'multinomial'));
            mw = fit_decoder(Wtr, ytr, nc, method);
            pw = apply_decoder(mw, Wte);
            for k = 1:numel(sq)
                res.(vn).win_pred((sq(k) - 1) * nc + yte(k)) = pw(k);
            end
            % permutation null (training window labels shuffled within sequence)
            if any(strcmp(vn, {'ica_ownz', 'multinomial'}))
                col = 1 + strcmp(vn, 'multinomial');
                g = findgroups(windows_seq(str, ltr, nc));
                for q = 1:opts.n_perm
                    rng(opts.seed + 1000 * q + f);
                    yp = ytr;
                    for gg = unique(g)'
                        ii = find(g == gg); yp(ii) = yp(ii(randperm(numel(ii))));
                    end
                    mp = fit_decoder(Wtr, yp, nc, method);
                    pp = apply_decoder(mp, Wte);
                    if isempty(perm_pred{q, col}); perm_pred{q, col} = zeros(n_seq * nc, 1); end
                    for k = 1:numel(sq); perm_pred{q, col}((sq(k) - 1) * nc + yte(k)) = pp(k); end
                end
            end
        end
    end
    for q = 1:opts.n_perm
        perm_acc(q, 1) = mean(perm_pred{q, 1} == win_true);
        perm_acc(q, 2) = mean(perm_pred{q, 2} == win_true);
    end
    labelled = lab_of_bin > 0;
    for v = 1:numel(variants)
        vn = variants{v};
        r = res.(vn);
        r.acc_bin = mean(r.bin_pred(labelled) == lab_of_bin(labelled));
        r.acc_win = mean(r.win_pred == win_true);
        r.conf_bin = confusionmat(lab_of_bin(labelled), r.bin_pred(labelled), 'Order', 1:nc) ./ accumarray(lab_of_bin(labelled), 1, [nc 1]);
        r.conf_win = confusionmat(win_true, r.win_pred, 'Order', 1:nc) ./ accumarray(win_true, 1, [nc 1]);
        res.(vn) = r;
        fprintf('  %-12s single 20 ms bins %.1f%% | 150 ms windows %.1f%%  (chance 20%%)\n', vn, 100 * r.acc_bin, 100 * r.acc_win);
    end
    fprintf('  assemblies per fold: %s; window null 95%%: ICA own z %.1f%%, multinomial %.1f%%\n', mat2str(n_asm'), ...
        100 * prctile(perm_acc(:, 1), 95), 100 * prctile(perm_acc(:, 2), 95));

    % all-data fit, for the figure: assembly profiles over classes and time
    mu = mean(Xb, 1); sd = std(Xb, 0, 1); sd(sd == 0) = 1;
    Z = (Xb - mu) ./ sd;
    A = assembly_templates(Z, 'ica', opts.seed);
    Rall = assembly_activation(Z, A.W);
    prof = zeros(size(A.W, 2), nc);
    for c = 1:nc; prof(:, c) = mean(Rall(lab_of_bin == c, :), 1)'; end
    tcourse = squeeze(mean(reshape(Rall, nb, n_seq, []), 2));     % bin x assembly
    if size(A.W, 2) == 1; tcourse = tcourse(:); end
    [~, pref] = max(prof, [], 2);
    pops(p) = struct('name', P.name, 'cluster_id', ids, 'n_assemblies_fold', n_asm, 'results', res, ...
        'perm_acc', perm_acc, 'W', A.W, 'lambda', A.lambda, 'lambda_max', A.lambda_max, 'profile', prof, ...
        'preferred_class', pref, 'timecourse', tcourse); %#ok<AGROW>
    fprintf('  all-data assemblies: %d; preferred class counts (A B C D grey): %s\n', size(A.W, 2), ...
        mat2str(accumarray(pref, 1, [nc 1])'));
end

results = struct('subject', subject, 'date', date, 'stimulus', stimulus_name, 'opts', opts, 'variants', {variants}, ...
    'variant_labels', {vlab}, 'classes', {cls}, 't', tc, 'bin_label', lab, 'fold', fold, 'trials', trials, 'pop', pops);
fig = plot_cv(results, tag);
if ~opts.save; return; end
base = sprintf('decoding_%s_assembly_cv', stimulus_name);
safe_save(fullfile(out_dir, [base '.mat']), struct('results', results), cfg);
save_vector_figures(fig, {base}, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function [Wf, y, sq] = windows(Fb, lb, sb, nc, use_sum)
% one row per (sequence, class) window: mean (or summed) features over the window's bins
key = (sb - 1) * nc + lb;
ok = lb > 0;
[u, ~, g] = unique(key(ok));
Fo = Fb(ok, :);
if use_sum
    Wf = splitapply(@(x) sum(x, 1), Fo, g);
else
    Wf = splitapply(@(x) mean(x, 1), Fo, g);
end
y = mod(u - 1, nc) + 1;
sq = floor((u - 1) / nc) + 1;
end

function s = windows_seq(sb, lb, nc)
key = (sb - 1) * nc + lb;
u = unique(key(lb > 0));
s = floor((u - 1) / nc) + 1;
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function f = plot_cv(r, tag)
np = numel(r.pop);
nc = numel(r.classes);
ccol = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655; 0.55 0.55 0.55];
cmap = interp1([0 0.6 1], [1 1 1; 0.165 0.471 0.839; 0.04 0.15 0.35], linspace(0, 1, 256));
f = figure('Color', 'w', 'Position', [30 30 1800 360 * np]);
tl = tiledlayout(f, np, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - assembly activation (PCA \\rightarrow Marchenko-Pastur \\rightarrow ICA) decodes A/B/C/D/grey on ' ...
    'held-out awake sequences? 10-fold CV, %d ms bins, LDA on activation vectors'], tag, round(1000 * r.opts.bin)), 'FontSize', 10.5);
for p = 1:np
    P = r.pop(p);
    pname = sprintf('%s (n = %d)', strrep(P.name, '_', ' '), numel(P.cluster_id));
    % 1. assembly profiles across classes
    ax = nexttile(tl);
    prof = P.profile;
    if ~isempty(prof)
        pz = (prof - mean(prof, 2)) ./ max(std(prof, 0, 2), eps);
        [~, o] = sortrows([P.preferred_class -max(prof, [], 2)]);
        imagesc(ax, pz(o, :)); colormap(ax, cmap); cb = colorbar(ax); cb.Label.String = 'activation (z across classes)';
        set(ax, 'XTick', 1:nc, 'XTickLabel', r.classes);
    end
    xlabel(ax, 'class'); ylabel(ax, 'assembly (sorted by preference)');
    title(ax, sprintf('%s: %d assemblies (all data)', pname, size(P.W, 2)), 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
    % 2. assembly time courses, coloured by preferred class
    ax = nexttile(tl); hold(ax, 'on');
    t = 1000 * r.t;
    for k = 1:size(P.timecourse, 2)
        plot(ax, t, P.timecourse(:, k), '-', 'Color', [ccol(P.preferred_class(k), :) 0.8], 'LineWidth', 1);
    end
    for e = 1:nc - 1
        x0 = 1000 * r.trials.element_start(e);
        xline(ax, x0, ':', r.classes{e}, 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    end
    xline(ax, 1000 * r.opts.grey, '-', 'Color', [0.7 0.7 0.7]);
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'mean activation strength');
    title(ax, 'assembly activation across the sequence (colour = preferred class)', 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
    % 3. accuracy per variant
    ax = nexttile(tl); hold(ax, 'on');
    acc = zeros(numel(r.variants), 2);
    for v = 1:numel(r.variants)
        acc(v, :) = 100 * [P.results.(r.variants{v}).acc_bin, P.results.(r.variants{v}).acc_win];
    end
    b = bar(ax, acc, 'BarWidth', 0.85, 'EdgeColor', 'none');
    b(1).FaceColor = [0.6 0.75 0.92]; b(2).FaceColor = [0.165 0.471 0.839];
    yline(ax, 20, ':', 'chance', 'FontSize', 7, 'HandleVisibility', 'off');
    plot(ax, b(2).XEndPoints(2), 100 * prctile(P.perm_acc(:, 1), 95), 'k_', 'MarkerSize', 12, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, b(2).XEndPoints(5), 100 * prctile(P.perm_acc(:, 2), 95), 'k_', 'MarkerSize', 12, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    for v = 1:numel(r.variants)
        text(ax, b(2).XEndPoints(v), acc(v, 2) + 2, sprintf('%.0f', acc(v, 2)), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    set(ax, 'XTick', 1:numel(r.variants), 'XTickLabel', r.variant_labels); xtickangle(ax, 30);
    ylim(ax, [0 105]); ylabel(ax, 'accuracy (%)');
    legend(ax, {'single 20 ms bins', '150 ms windows'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    title(ax, 'held-out accuracy (- = window null 95%)', 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
    % 4-5. confusion (ICA, own z): windows and single bins
    tt = {'conf_win', '150 ms windows'; 'conf_bin', 'single 20 ms bins'};
    for q = 1:2
        ax = nexttile(tl);
        C = 100 * P.results.ica_ownz.(tt{q, 1});
        imagesc(ax, C, [0 100]); colormap(ax, cmap); axis(ax, 'square');
        for i = 1:nc
            for j = 1:nc
                text(ax, j, i, sprintf('%.0f', C(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [1 1 1] * (C(i, j) > 55));
            end
        end
        set(ax, 'XTick', 1:nc, 'XTickLabel', r.classes, 'YTick', 1:nc, 'YTickLabel', r.classes);
        xlabel(ax, 'decoded'); ylabel(ax, 'true');
        title(ax, sprintf('ICA, own z: %s', tt{q, 2}), 'FontWeight', 'normal', 'FontSize', 8.5);
        style(ax);
    end
end
end
