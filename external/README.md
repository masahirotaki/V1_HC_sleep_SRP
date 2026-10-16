# external/VR_NPX_analysis

Unmodified copies of the functions this repository needs from
`C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis` (copied 2026-10-13),
kept in their original sub-folder layout (includes buzcode / FMAToolbox helpers and their mex files).

Entry points used here:
- `behavioural_state_detector/detect_behavioural_states_masa.m` — sleep staging
- `LFP_analysis/FindRipples_masa.m` — ripples
- `LFP_analysis/DetectSpindles_masa.m` — spindles
- `LFP_analysis/DetectSlowWaves_masa.m` — slow oscillations / UP-DOWN states
- `utility/select_clusters.m`, `utility/create_cluster_selection_params.m` — unit quality
- `spikes/analysis/psthAndBA.m` — PSTHs (for later analyses)
- `utility/save_all_figures.m` — figure saving (.fig + .pdf, or .svg + raster .pdf), called via `utils/save_figures.m`

Change behaviour through the wrappers in `sleep/` and `config/project_config.m` rather than editing these copies.
