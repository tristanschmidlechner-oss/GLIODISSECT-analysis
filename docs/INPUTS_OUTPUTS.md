# Inputs and outputs

## Normative arrays

`Y` is voxels × controls, real single/double; NaN denotes missing observations.
Inf is rejected, including overflow when converting to single. Single precision
loading followed by double precision sufficient-statistic OLS matches the validated
G-model workflow. All columns and G values must share control order and all rows
must share spatial geometry. No patient belongs in the training cohort.

`anatomy` is logical and the same size as Y. Study G-model anatomy was individual
registered GM+WM >=0.8, sampled on the response grid with SPM interpolation 1,
both tissue values finite. This is not a cohort-intersection mask. `G` is a finite
control vector: median finite native-T1 oxygen-normalized step-1 response in the
established individual G mask. `pone_g_scalar` calculates that median from a
caller-supplied mask; it does not construct or validate that mask. Patient G-mask
tumor exclusions are caller responsibilities. No cohort centering/scaling is used.

`pone_fit_conventional(Y,valid,cfg)` instead takes each map's explicit logical
support. The conventional reference uses finite-map coverage, without
adding the G-model tissue-coverage rule. Both response references use %BOLD/mmHg
when applied to oxygen-normalized response maps. Timing arrays are in volumes
(multiply by 1.8 for seconds); raw response is %BOLD. G must use the modeled
oxygen-normalized step-1 response units. z, R² and Q² are dimensionless.

Default counts are controls=20, full_required=16, training_required=16,
q2_required=17. For another cohort supply all four fields, as shown in
[example_configurations](../examples/example_configurations.m). No coverage
rounding or automatic tuning is performed. At least four controls and three
training observations are required for residual degrees of freedom.

## Patient application

`pone_apply_g(observed,G,mask,model)` and
`pone_apply_conventional(observed,mask,model)` require all arrays to have identical
shape and already occupy the patient's grid. They enforce different `kind` tags.
No implicit warp or flattening is performed. Models may be used on their original
grid, or the caller must resample the relevant fields and logical support explicitly.
For the established native workflow, continuous reference maps use fourth-order
interpolation and support nearest-neighbour. This package does not perform that warp.

G-adjusted raw residual SD is floored after resampling. Conventional sample SD is
floored before resampling and must remain finite and strictly positive afterward;
small positive values are retained without clipping/reflooring. Therefore a native
conventional SD can be below the original floor and extreme z can occur. Outputs:
expected, difference, z, signed_abnormal (-1/0/+1, strict ±2), valid_mask and a
nonpositive-SD count. Invalid deviations and signed classes are NaN.

## LOOCV outputs

`pone_write_loocv` requires a new directory. `LOOCV_RESULTS.csv` contains one row
of cross-validated aggregates; `NORMATIVE_MODEL_RESULTS.csv` one in-sample row.
`LOOCV_subjects.csv` contains ordinal fold metrics, `LOOCV_coverage.csv` validity
counts, and `summary.json` separate statistics plus the exact configuration.
`voxelwise_metrics.mat` contains both Q² maps, PRESS, held-out TSS/counts,
eligibility, and explicitly named in-sample R². MAT is v7.3. No input arrays are
serialized. Nevertheless real-data voxelwise/fold outputs are private and must
not be added to the distributable ZIP. Zero eligible scores yield NaN metrics;
no positive training residual support raises an error.

## Temporal/smoothing inputs

`pone_thresholds(input,mask,pair)` takes `bBOLD` (X×Y×Z×160 filtered data),
`brainLag` (initial lag in volumes), `durationRA` (five phase durations in volumes),
`initialGlobalLag`, `rowScans=1:160`, and `shiftedO2` (aligned mmHg oxygen vector
or callback accepting the refined global lag). See the synthetic input factory.
The logical mask must share the spatial shape. Acquisition: TR=1.8 s, two 60-s
steps separated by 40 s. Other acquisitions need adaptation and validation.

`pone_smoothing_workflow(paths,mask,endTidalFile,eventsFile,cfg,newOutputDir)` takes
160 explicitly ordered, aligned 3-D NIfTI paths, and the external release's
RespirAct 23-field EndTidal text/Events spreadsheet formats. Oxygen is field 6,
timestamps field 1 in milliseconds; event times are column 1 and scan indices
column 2. See the external release's input guide for complete acquisition requirements.
Spatial 4/8-mm alternatives require unsmoothed data; temporal 8/16 alternatives
use already 6-mm-smoothed data (`input_is_unsmoothed=false`). Re-smoothing an
already smoothed input is not equivalent. A 6-mm,12-volume baseline can be run
with either appropriately supplied input kind. Input zeros become NaN before
filtering, as in the directory workflow. Mask construction and spatial
preprocessing are external responsibilities; do not relabel a supplied mask as G anatomy.

The driver saves sensitivity.mat, convergence.csv and 13 released-map NIfTIs.
It preserves intermediate smoothed data under its new output directory.
Partial outputs remain after failure and are never overwritten on retry.
