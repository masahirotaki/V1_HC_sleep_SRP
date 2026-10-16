Copied unchanged on 2026-10-13 from
`V1_MEC_acute_MAT\dim_high_functions\Toolboxes_customised\DataHigh-master` (the customised DataHigh
used by `V1_MEC_acute_MAT\dimension_reduction\FA_continuous_test.m`), so the FA here matches that code.
DataHigh: Cowley, Kaufman, Butler, Churchland, Ryu, Shenoy & Yu (2013), GPL v2 (COPYING.txt).

- `reducedims_EH.m`  (root version, returns FA `params`); alg 1 = PCA, 3 = FA. Only PCA/FA are used here;
  the PPCA/LDA/GPFA branches need DataHigh files that were not copied.
- `fastfa.m`, `fastfa_estep.m`, `cosmoother_fa.m`  (gpfa/core_twostage)
- `orthogonalize.m`, `assignopts.m`, `logdet.m`  (gpfa/util)
