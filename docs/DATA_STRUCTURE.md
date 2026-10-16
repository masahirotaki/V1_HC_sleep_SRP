# Ellie's data (Y:\Ellie\DATA)

Read-only, except for our outputs: they are written into each recording's own analysis folder,
`SUBJECTS\<subject>\analysis\<date>\<recording>\` (and its `detection\` and `Decoding\` subfolders). Existing files
there are never overwritten, except `behavioural_state.mat`, which our pipeline replaces after
backing up Ellie's version as `behavioural_state_Ellie.mat`
(`check_output_path`, `check_overwrite`, `cfg.output.own_files`).

## Session table
`Y:\Ellie\DATA\Ellie_session_stimuli_table.xlsx`, one sheet per animal:
`animal, date, gs_number, stimulus_type, probe_number, probe_hemisphere, session_id, ...`
`session_id` → ephys sub-folder, `gs_number` → SpikeGLX g-index. Resolved by `get_session_paths`.

## Folders (per subject)
```
SUBJECTS/<subj>/
  ephys/<date>/<date>_<sid>/<subj>_<date>_<sid>_g<gs>/<..>_imec0/
      *_tcat.imec0.lf.bin/.meta   CatGT LFP: 2.5 kHz, 385 ch, global CAR, 0-500 Hz
      *_t0.imec0.ap.bin/.meta     raw AP (30 kHz)
      *_ap_kilosortChanMap.mat    xcoords / ycoords (depth, um) / kcoords (shank, 1-based)
  ephys/<date>/probe0/            SpikeInterface sorting of all recordings of the day (concatenated)
  stimuli/<date>/<sid>_g<gs>/     Bonsai csv (video mobility, sync)
  analysis/<date>/
      depths_from_PSD.mat         .shank<k>.CA1_depth_PSD / L5_depth_PSD / surface_depth_PSD  (NPX1: no shank level)
      earliest_V1sink_CSD.mat     struct array .shank_id / .best_channel_this_shank / .best_depth_this_shank (L4)
      <recording>/                e.g. Sleep_Box_1, GAVNIK_ABCD_1 (repeated names -> <name>_<k>)
          extracted_behaviour.mat     Behaviour
          extracted_clusters_ks4.mat  clusters_ks4
          extracted_task_info.mat     Task_info (visual stimulus recordings: stim_onset, stim_orientation)
          ripples.mat, spindles.mat, slowwaves.mat  (Ellie's earlier sleep outputs)
          behavioural_state.mat, extracted_LFP.mat, behavioural_state/sleep_stages .fig/.pdf  (ours)
          detection/                  extracted_*_events*.mat, detection_summary.csv, figures (ours)
          Decoding/                   decoding_* results / figures (ours, decoding/run_abcd_decoding)
```

## Visual stimulus timing (`extracted_task_info.mat` → `Task_info`, GAVNIK ABCD)
- `stim_onset` = one entry per **element** onset (SpikeGLX s), not per sequence; `stim_orientation` per element.
  M00069 20250929 GAVNIK_ABCD_1: 800 onsets = 200 sequences × A,B,C,D = 10°,100°,40°,130°, always in that order.
- Elements are back to back (SOA ~150 ms); no offsets are stored. D offset = D onset + element duration;
  then ~1.52 s grey until the next A. Pauses of ~30 s split the 200 sequences into 4 blocks of 50.
- `pd_on/pd_off.sglxTime` (400 each): the photodiode flips at every element (on = A, C; off = B, D).
- Parsed by `get_abcd_trials`.

## Time base
All of these share SpikeGLX (imec probe) time in seconds, 0 = first sample of the recording:
LFP sample index / imSampRate, `Behaviour.sglxTime`, `clusters_ks4.spike_times`.

## Behaviour (`extracted_behaviour.mat` → `Behaviour`)
60 Hz, aligned to SpikeGLX through the Bonsai–NIDQ sync pulse.
- `tvec` = `sglxTime` (s); starts a few seconds before / ends after the ephys file
- `mobility` — video pixel change (a.u.); `X`, `Y` — tracked position (px)

Immobility (`load_session_behaviour`, Ellie's method): `|diff(movmean(mobility, 1 s))| > 2000` = moving,
pauses < 1 s filled in. `mobility_zscore` = z-scored version, used by the ripple/spindle detectors (exclude z > 3).

## Clusters (`extracted_clusters_ks4.mat` → `clusters_ks4`)
- `cluster_id` (1-based), `spike_times` (s), `spike_id` (cluster of each spike)
- `peak_channel` → channel-map row → `ycoords` depth, `kcoords` shank
- quality metrics; Ellie's selection: `create_cluster_selection_params('sorting_option','ellie')`
- region (`load_session_clusters`): V1 = L5 depth −330…+700 um, CA1 = CA1 depth ±150 um, per shank

## Shanks
NPX2.0 4-shank (NP2013) animals record two shanks: {2,3} (M00069, M00088, M00113) or {1,4}
(M00087, M00096, M00100). Ripples use shank 4, else 3 (most anterior CA1). M00013/M00014 are NPX1.0.
