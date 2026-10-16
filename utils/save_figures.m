function save_figures(figs, names, folder, cfg)
% SAVE_FIGURES  Save each figure with save_all_figures (.fig + .pdf) into FOLDER.
%
% Figures whose scatter / marker-only plots hold more than cfg.plot.svg_min_points
% points use save_all_figures' SVG option (.svg + rasterised .pdf) to keep files small.
% save_all_figures saves every open figure and names files by the figure Name, so
% each figure is saved on its own with all other figures hidden from it. Figures are
% closed after saving.

check_output_path(folder, cfg);
% save_all_figures silently redirects into <folder>\figures\png_figs if that exists
assert(~exist(fullfile(folder, 'figures'), 'dir'), ...
    'save_figures: %s contains a figures\\ subfolder; save_all_figures would write there instead', folder);
if ~exist(folder, 'dir'); mkdir(folder); end
for i = 1:numel(figs)
    fig = figs(i);
    name = regexprep(names{i}, '[^\w\-]', '_');
    use_svg = n_marker_points(fig) > cfg.plot.svg_min_points;
    outputs = {'.fig', '.pdf'};
    if use_svg; outputs{end + 1} = '.svg'; end %#ok<AGROW>
    for k = 1:numel(outputs)
        check_overwrite(fullfile(folder, [name outputs{k}]), cfg);
    end

    all_figs = findobj(groot, 'Type', 'figure');
    others = all_figs(all_figs ~= fig);
    set(others, 'HandleVisibility', 'off');
    restore = onCleanup(@() set(others(isvalid(others)), 'HandleVisibility', 'on'));
    set(fig, 'Name', name);
    set(fig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); % saved .fig opens visible
    set(groot, 'CurrentFigure', fig);
    if use_svg
        save_all_figures(folder, [], 'SVG_option', true);
    else
        save_all_figures(folder, []);
    end
    if isvalid(fig); close(fig); end
    clear restore

    missing = outputs(~cellfun(@(ext) isfile(fullfile(folder, [name ext])), outputs));
    if ~isempty(missing)
        warning('save_figures:missing', 'Not written: %s (%s) in %s', name, strjoin(missing, ', '), folder);
    end
end
end


function n = n_marker_points(fig)
n = sum(arrayfun(@(h) numel(h.XData), findobj(fig, 'Type', 'scatter')));
markers = findobj(fig, 'Type', 'line', 'LineStyle', 'none');
n = n + sum(arrayfun(@(h) numel(h.XData), markers));
end
