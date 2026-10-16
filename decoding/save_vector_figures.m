function save_vector_figures(figs, names, folder, cfg)
% SAVE_VECTOR_FIGURES  Save each figure as .fig + vector .pdf (save_figures), never as a
% raster: dense scatter plots are not switched to SVG + rasterised PDF, and no PNG is
% written. Large raster/scatter figures make bigger, slower PDFs.

cfg.plot.svg_min_points = inf;
save_figures(figs, names, folder, cfg);
end
