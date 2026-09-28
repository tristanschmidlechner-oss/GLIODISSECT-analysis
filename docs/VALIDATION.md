# Validation

## Synthetic tests

`run_pone_tests` exercises independent backslash-OLS and Q² oracles, training/held-out
isolation, coverage boundaries, residual floors, pooled weighting, distinct patient
models, configurable cohort sizes, invalid-input rejection and CSV/JSON/MAT output.
Adapter tests compare all returned GLIODISSECT fields exactly across 12 synthetic
ensembles, including invalid voxels (396 field comparisons). Six finite-data
ensembles also compare common fields with `pone_study_core` (216 comparisons).
SPM tests compare 4/6/8-mm smoothing directly with SPM, exercise 8/12/16-volume
filtering, and round-trip the GLIODISSECT NIfTI output. The complete smoothing
workflow is tested in five configurations with synthetic NIfTI, EndTidal and Events
inputs; all release fields are compared exactly to a direct release pipeline.
The LOOCV implementation also matched the study driver exactly on synthetic data.

The [recorded MATLAB test results](../validation/test_status.json) document
synthetic validation in MATLAB R2025b on 28 September 2026.
Run `run_pone_tests` to verify the package in your MATLAB installation.
Examples use fixed seeds; their numbers are not study targets. Synthetic tests
do not establish biological accuracy or clinical validity.

## Study core and missing-data conventions

Study reference aggregates were generated with the study-analysis workflow.
The included private `pone_study_core` implements its temporal-core equations
and support conventions for comparison with the analysis adapter. It does not
exclude nonfinite BOLD series before landmark estimation and returns DTP2 without
an output mask. GLIODISSECT v1.0.0 and the adapter exclude invalid input series
from output maps and lag references and mask DTP2. The finite-data comparisons
above verify exact agreement of common output fields on the tested synthetic
inputs; they do not establish equivalence for missing data or an entire cohort.

## Study reference results

Study LOOCV validation in MATLAB R2025b/SPM12 reproduced 17 aggregate metrics
and all 20 fold rows. The public reference JSON contains aggregate results only.
Threshold/convergence validation included exact baseline comparisons in all
20 controls. The sensitivity reference tables summarize the completed study
analyses; they are distinct from the synthetic test outputs.

No participant processing was repeated to package this supplement. Exact study
reproduction requires private imaging inputs, masks, G scalars and registration
provenance. Registration and regional atlas/tissue/G sensitivity workflows are
outside the implemented scope.
