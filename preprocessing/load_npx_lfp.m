function LFP = load_npx_lfp(lf_bin_file, bin_channels, cfg)
% LOAD_NPX_LFP  Read selected channels from a SpikeGLX lf.bin (read-only),
% convert to uV and downsample to cfg.lfp.target_fs.
%
%   LFP = load_npx_lfp(paths.lf_bin_file, layout.chanMap(rows), cfg)
%
% bin_channels  1-based channel indices in the binary file (layout.chanMap(row)).
% LFP.data      [nChannels x nSamples] double, uV
% LFP.tvec      [1 x nSamples] seconds from the first sample of the file, i.e.
%               SpikeGLX/imec time - the same clock as Behaviour.sglxTime and
%               cluster spike_times.
%
% CatGT "tcat" LFP files are already global-CAR'd and low-passed (500 Hz), so
% cfg.lfp.car = 'auto' only re-references raw (non-tcat) files, using the median
% across all ephys channels. Before decimation a zero-phase low-pass at 80% of the
% new Nyquist is applied to the whole trace (no chunk-edge artefacts).

meta = read_spikeglx_meta(regexprep(lf_bin_file, '\.bin$', '.meta'));
n_saved = str2double(meta.nSavedChans);
fs_in = str2double(meta.imSampRate);
ap_lf_sy = str2double(strsplit(meta.snsApLfSy, ','));
n_ephys = sum(ap_lf_sy(1:2));
file_info = dir(lf_bin_file);
n_samp = floor(file_info.bytes / (2 * n_saved));
uv_per_bit = microvolts_per_bit(meta);

switch cfg.lfp.car
    case 'auto'; do_car = ~contains(lf_bin_file, 'tcat');
    otherwise;   do_car = logical(cfg.lfp.car);
end

factor = max(1, round(fs_in / cfg.lfp.target_fs));
chunk = round(cfg.lfp.chunk_sec * fs_in);
n_chunks = ceil(n_samp / chunk);
raw = zeros(numel(bin_channels), n_samp, 'single');

fprintf('Loading %d LFP channels from %s (%.1f min, %d chunks, CAR=%d)\n', ...
    numel(bin_channels), lf_bin_file, n_samp / fs_in / 60, n_chunks, do_car);
fid = fopen(lf_bin_file, 'r');
closer = onCleanup(@() fclose(fid));
pos = 0;
for c = 1:n_chunks
    n = min(chunk, n_samp - pos);
    block = fread(fid, [n_saved n], '*int16');
    if do_car
        x = single(block(1:n_ephys, :));
        x = x - median(x, 1);
        raw(:, pos + (1:n)) = x(bin_channels, :);
    else
        raw(:, pos + (1:n)) = single(block(bin_channels, :));
    end
    pos = pos + n;
    if mod(c, 10) == 0 || c == n_chunks; fprintf('  chunk %d/%d\n', c, n_chunks); end
end

fs_out = fs_in / factor;
keep = 1:factor:n_samp;
LFP.data = zeros(numel(bin_channels), numel(keep));
if factor > 1
    [z, p, k] = butter(8, 0.8 * (fs_out / 2) / (fs_in / 2), 'low');
    [sos, g] = zp2sos(z, p, k);
end
for ch = 1:numel(bin_channels)
    x = double(raw(ch, :)) * uv_per_bit;
    if factor > 1; x = filtfilt(sos, g, x); end
    LFP.data(ch, :) = x(keep);
end

LFP.tvec = (keep - 1) / fs_in;
LFP.fs = fs_out;
LFP.bin_channels = bin_channels(:)';
LFP.source_file = lf_bin_file;
LFP.car_applied = do_car;
LFP.uv_per_bit = uv_per_bit;
end


function c = microvolts_per_bit(meta)
range_max = str2double(meta.imAiRangeMax);
if isfield(meta, 'imMaxInt'); max_int = str2double(meta.imMaxInt); else; max_int = 512; end
if isfield(meta, 'imChan0lfGain')
    gain = str2double(meta.imChan0lfGain);
elseif isfield(meta, 'imChan0apGain') % NPX2: LFP is derived from the AP band
    gain = str2double(meta.imChan0apGain);
else % fall back to the imro table
    entries = regexp(meta.imroTbl, '\(([^\)]*)\)', 'tokens');
    first = str2double(strsplit(entries{2}{1}, ' '));
    if numel(first) == 6; gain = first(5); else; gain = 80; end % NPX1 lfGain / NPX2 fixed gain
end
c = range_max / max_int / gain * 1e6;
end
