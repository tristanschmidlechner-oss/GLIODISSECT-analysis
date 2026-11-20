# Validation

## Synthetic tests

`run_pone_tests` exercises independent backslash-OLS and Q² oracles, training/held-out
isolation, coverage boundaries, residual floors, pooled weighting, distinct patient
models, configurable cohort sizes, invalid-input rejection and CSV/JSON/MAT output.
Adapter tests compare all returned GLIODISSECT fields exactly across 12 synthetic
ensembles, including invalid voxels (396 field comparisons). Three additional
oxygen cases cover zero and reversed oxygen differences and the oxygen callback
interface (99 exact release-field comparisons).
The four supplementary whole-stimulus-average fields are checked against an
independent phase-weight calculation on all 15 inputs (60 field checks, absolute
tolerance 1e-11). That oracle uses released lag, support and oxygen-difference
outputs, without calling the supplement's averaging helper. It checks missing
output support separately from finite arithmetic, including zero oxygen division.
SPM tests compare 4/6/8-mm smoothing directly with SPM, exercise 8/12/16-volume
filtering, and round-trip the GLIODISSECT NIfTI output. The complete smoothing
workflow is tested in five configurations with synthetic NIfTI, EndTidal and Events
inputs; all release fields are compared exactly to a direct release pipeline.

The [MATLAB test results](../validation/test_status.json) record the synthetic
checks, software versions, release revision and standalone example run.
Run `run_pone_tests` to verify the package in your MATLAB installation.
Examples use fixed seeds; their numbers are not study targets. Synthetic tests
do not establish biological accuracy or clinical validity.

## Output support and interpretation

GLIODISSECT v1.0.0 and the adapter exclude nonfinite BOLD series and nonfinite
initial lags from eligible input support. Lag references and output maps require
a finite refined lag within that support; DTP2 is masked like the other maps.
Nonfinite map outputs are represented as NaN. The returned analysis mask and
validity mask have distinct meanings: requested spatial support and eligible
output support, respectively.

Exact release-field agreement applies to the tested synthetic inputs, including
the specified invalid-input cases. It does not establish equivalence for every
possible missing-data pattern or reproduce an entire cohort. The supplement's
four whole-stimulus-average fields extend the released output set; their
baseline and oxygen-normalization definitions differ from steady-state responses
as described in [sensitivity and convergence](SENSITIVITY.md).

## Study reference results

The final 20-control rerun in MATLAB R2025b/SPM12 reproduced all 19 numeric
aggregate values in `reference/study_loocv_aggregates.json` (14 cross-validated
and 5 in-sample values) and all 20 fold rows (14 shared numeric fields per row).
The aggregate values include counts as well as performance measures. Maximum
absolute aggregate error was 1.4210854715202004e-14, within the comparison
tolerance of 1e-10 * max(1, abs(reference)). This rerun refit the models from
the recorded participant MNI response maps, anatomy eligibility and G scalars;
it did not recompute registration or G-mask derivation. The public reference
JSON contains aggregate results only.
Threshold/convergence validation included exact baseline comparisons in all
20 controls. The sensitivity reference tables summarize the completed study
analyses; they are distinct from the synthetic test outputs.

Exact study reproduction requires private imaging inputs, masks, G scalars and registration
provenance. Registration and regional atlas/tissue/G sensitivity workflows are
outside the implemented scope.
