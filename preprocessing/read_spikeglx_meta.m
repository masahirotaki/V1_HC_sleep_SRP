function meta = read_spikeglx_meta(meta_file)
% READ_SPIKEGLX_META  Parse a SpikeGLX .meta file into a struct of strings.
% Keys starting with '~' (e.g. ~imroTbl) are stored without the '~'.

txt = fileread(meta_file);
lines = regexp(txt, '\r?\n', 'split');
meta = struct();
for i = 1:numel(lines)
    kv = regexp(lines{i}, '^~?([^=]+)=(.*)$', 'tokens', 'once');
    if isempty(kv); continue; end
    meta.(matlab.lang.makeValidName(strtrim(kv{1}))) = strtrim(kv{2});
end
end
