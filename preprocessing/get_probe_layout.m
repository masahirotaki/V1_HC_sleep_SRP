function layout = get_probe_layout(paths, cfg)
% GET_PROBE_LAYOUT  Channel geometry plus CA1 / L5 / L4 channels on every shank.
%
% depths_from_PSD.mat      NPX2: .shank<k>.CA1_depth_PSD / L5_depth_PSD / surface_depth_PSD
%                          NPX1: the same fields without the shank level
% earliest_V1sink_CSD.mat  NPX2: struct array with .shank_id / .best_depth_this_shank
%                          NPX1: .overall_best_halfmax_depth
% Channels are returned as channel-map row indices (1-based); the matching
% column of the binary file is layout.chanMap(row). Cluster peak_channel uses
% the same row indexing.

cm = load(paths.chanmap_file);
layout.xcoords = cm.xcoords(:);
layout.ycoords = cm.ycoords(:);
layout.kcoords = cm.kcoords(:);
layout.chanMap = cm.chanMap(:);

psd = load(paths.psd_file);
psd = psd.depths_from_PSD;
if exist(paths.csd_file, 'file')
    csd = load(paths.csd_file);
    csd = csd.earliest_V1sink_CSD;
else
    warning('No earliest_V1sink_CSD.mat for %s %s: L4 depth unavailable.', paths.subject, paths.date);
    csd = [];
end

fn = fieldnames(psd);
shank_fields = fn(startsWith(fn, 'shank'));
if isempty(shank_fields) % NPX1.0: single shank
    shank_ids = unique(layout.kcoords)';
    shank_psd = {psd};
else
    shank_ids = cellfun(@(f) str2double(erase(f, 'shank')), shank_fields)';
    shank_psd = cellfun(@(f) psd.(f), shank_fields, 'UniformOutput', false);
end

for i = 1:numel(shank_ids)
    s.id            = shank_ids(i);
    s.CA1_depth     = shank_psd{i}.CA1_depth_PSD;
    s.L5_depth      = shank_psd{i}.L5_depth_PSD;
    s.surface_depth = shank_psd{i}.surface_depth_PSD;
    s.L4_depth      = csd_depth(csd, s.id);
    s.CA1_channel   = channel_at_depth(layout, s.id, s.CA1_depth);
    s.L5_channel    = channel_at_depth(layout, s.id, s.L5_depth);
    s.L4_channel    = channel_at_depth(layout, s.id, s.L4_depth);
    s.V1_range      = s.L5_depth  + cfg.depth_ranges.V1;
    s.CA1_range     = s.CA1_depth + cfg.depth_ranges.CA1;
    s.L4_range      = s.L4_depth  + cfg.depth_ranges.L4;
    layout.shank(i) = s;
end

%% Detection channels
ids = [layout.shank.id];
pick = find(ismember(cfg.ca1_shank_preference, ids), 1);
if isempty(pick)
    layout.ripple_shank = ids(end);
    warning('None of shanks %s recorded; using shank %d for CA1.', mat2str(cfg.ca1_shank_preference), ids(end));
else
    layout.ripple_shank = cfg.ca1_shank_preference(pick);
end
if isempty(cfg.l5_shank)
    layout.l5_shank = layout.ripple_shank;
else
    layout.l5_shank = cfg.l5_shank;
end
layout.ripple_channel = layout.shank(ids == layout.ripple_shank).CA1_channel;
layout.l5_channel     = layout.shank(ids == layout.l5_shank).L5_channel;

fprintf('Shanks with depth info: %s\n', mat2str(ids));
fprintf('Ripple channel: shank %d, row %d (depth %g um)\n', layout.ripple_shank, layout.ripple_channel, layout.ycoords(layout.ripple_channel));
fprintf('L5 channel:     shank %d, row %d (depth %g um)\n', layout.l5_shank, layout.l5_channel, layout.ycoords(layout.l5_channel));
end


function depth = csd_depth(csd, shank_id)
depth = NaN;
if isempty(csd); return; end
if isfield(csd, 'shank_id')
    k = find([csd.shank_id] == shank_id, 1);
    if ~isempty(k); depth = csd(k).best_depth_this_shank; end
elseif isfield(csd, 'overall_best_halfmax_depth')
    depth = csd(1).overall_best_halfmax_depth;
end
end


function row = channel_at_depth(layout, shank_id, depth)
% Nearest channel to DEPTH on the shank (first column if two share the depth).
row = NaN;
if isnan(depth); return; end
on_shank = find(layout.kcoords == shank_id);
[~, k] = min(abs(layout.ycoords(on_shank) - depth));
row = on_shank(k);
end
