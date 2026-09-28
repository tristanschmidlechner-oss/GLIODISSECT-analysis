# Sensitivity and convergence

The supplement runs bounded DTP/DTR estimation with
10/90, 5/95 or 20/80 threshold pairs. The release does not expose threshold or
iteration-audit arguments. The adapter combines its public-core validity and
output conventions with separately named, private instrumented estimators.
It preserves author notices, indices, fixed windows, fallbacks, convergence
criteria and iteration bounds. No release file is modified. Baseline equivalence
is tested against GLIODISSECT v1.0.0 and the study core on synthetic data.

The returned `audit` has four records in DTP1,DTR1,DTP2,DTR2 order. Each record's
uint16 `rows` is voxels × 8 in MATLAB linear order: attempts, completed refinements,
criterion satisfied, bounded-unsatisfied flag, fallback bitmask, invalid crossing
or nonfinite criterion flag, caught error, skipped input/window. DTP may complete
up to 17 refinements and DTR up to 15. The summarizer uses the completed count to
identify reaching the bound; it does not confuse that with the bounded-unsatisfied flag.

DTP criterion: both successive baseline and step ratios finite and in
[0.9999,1.0001]. DTR: absolute baseline ratio minus one strictly <0.0001.
Evaluable counts require finite output and exclude caught/skipped estimates.
Unsatisfied percentages use evaluable finite estimates; fallback rates use mask
voxels. NaN/Inf/finite counts are reported separately and sum to mask size.
The final-index validity check uses the supported 160-volume acquisition.

Fallback bits (zero-based) record estimator events:

| Bit | DTP | DTR |
|---|---|---|
| 0 | Initial high crossing fallback | Initial low return crossing fallback |
| 1 | Initial low crossing fallback | Refined low return crossing fallback |
| 2 | Refined high crossing fallback | Recovery baseline-start window rule |
| 3 | Refined low crossing fallback | Final recovery-start window rule |
| 4 | Inverted step-window rule | No high return crossing retained |

DTP crossing mask=15/window mask=16; DTR crossing mask=19/window mask=12.
Flags are encounters during estimation, not proof that the final index uses a
fallback. The DTR high-crossing bit is cleared when a crossing is found in that
initial search. Numerical criterion satisfaction is not physiological equilibration.

`pone_compare_maps(...,'threshold')` uses the study field-specific
zero-to-NaN export convention and float32 conversion for threshold summaries.
Valid zero DTP2/relative-arrival/refined-lag values remain. This study-summary convention
is explicit and differs from the release's general finite-zero preservation.
The adapter's maps themselves follow the release; no zero removal occurs there.
`'smoothing'` accepts already exported maps and their explicit validity. For exact
study smoothing summaries the caller must supply the study exported support.
Threshold magnitude change uses medians of absolute responses; smoothing uses
absolute medians, according to their respective summary definitions. Timing remains volumes.
Threshold summaries include all 14 timing/steady-state/whole-average
quantities; smoothing summaries include seven timing/normalized-response quantities.
Whole-average responses use
the three-segment baseline and unguarded oxygen division; this differs
from the shared first-step steady-state baseline and its oxygen positivity rule.
Nonfinite outputs are invalid. The four whole-average maps remain in MAT outputs;
the external release NIfTI writer writes its own 13-map set.

Smoothing follows SPM spatial smoothing to float32, then zero-to-NaN, linear
detrending with restored mean, the release's 0.125-Hz smooth-box filter, robust
local smoothing, oxygen alignment and lag estimation. Spatial alternatives keep
12 volumes; temporal alternatives keep 6 mm. Joint alternatives are rejected.
The driver uses GLIODISSECT v1.0.0 for preprocessing and initial lag estimation.

The study reference tables summarize 40 threshold runs, 20 baseline audits and
80 smoothing runs. They are aggregate results, not outputs of the synthetic tests.
The threshold reference table also contains ICC(2,1) across paired subject
medians for baseline versus alternative thresholds. This quantifies threshold
sensitivity; it is not a comparison of steady-state and average maps. These
cohort ICC summaries are reference values; the executable summary function
returns per-subject map comparisons.

Regional atlas-based smoothing, tissue/G sensitivity, registration and regional
recurrence maps require separate workflows. This package provides numerical
outputs for downstream analysis; manuscript figures are not generated.
