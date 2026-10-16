function save_fig_formats(fig, out_dir, name)
% SAVE_FIG_FORMATS  Save a figure as .png (preview), vector .pdf and MATLAB .fig, then close it.
%
% The .fig opens visible even though figures are drawn with 'Visible','off'.
% The PDF uses exportgraphics ContentType 'vector', so lines, text and shading
% stay editable (e.g. in Illustrator).

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
base = fullfile(out_dir, name);
exportgraphics(fig, [base '.png'], 'Resolution', 130);
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
set(fig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
savefig(fig, [base '.fig'], 'compact');
close(fig);
end
