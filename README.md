# GLIODISSECT analysis supplement

Public manuscript-analysis supplement for PONE-D-26-25851, accompanying
GLIODISSECT v1.0.0. It provides executable analysis code and fully synthetic
examples. Exact study numbers require private study inputs, which are not included.

## Installation and synthetic quick start

Use MATLAB R2025b, Signal Processing Toolbox, Statistics and Machine Learning
Toolbox, SPM12, and the separate [GLIODISSECT v1.0.0 release](https://github.com/tristanschmidlechner-oss/GLIODISSECT/releases/tag/v1.0.0)
(commit `72fe3a4b511388f233dbd79b9afda456b17f3882`). No dependency is bundled.
Edit the three installation paths below. Do not add the supplement recursively.

```matlab
restoredefaultpath
supplementDir = '/path/to/GLIODISSECT-analysis';
gliodissectDir = '/path/to/GLIODISSECT';
spmDir = '/path/to/spm12';
addpath(supplementDir)
setup_pone_analysis(gliodissectDir,spmDir)
results = run_pone_tests();
out = fullfile(tempdir,'pone_synthetic_example'); % must not exist
R = example_synthetic(out);
```

The example writes synthetic LOOCV CSV/JSON/MAT outputs, both reference-model
definitions and their deviations, three threshold results, convergence counts
and threshold comparisons. Tests also exercise spatial smoothing with SPM and
temporal smoothing spans and the complete smoothing driver. To retain its synthetic
files, run `example_synthetic_smoothing(anotherNewOutputDir)`. Output files made from users' data remain private.

## Scope

| Component | Implemented behavior | Limits |
|---|---|---|
| G-adjusted reference | Voxelwise OLS, residual SD, fixed positive-SD P5 floor, patient application | Prepared aligned responses, individual anatomy and a separately calculated response-derived G scalar required |
| Conventional reference | Mean/sample SD, map-specific support, P5 floor before resampling, patient application | Distinct reference; caller supplies aligned/warped arrays |
| LOOCV | 20 folds of 19 controls; training-only masks/floors; matched mean comparator; Q² and subject/pooled RMSE | Other cohort sizes require explicit integer coverage counts; no external validation |
| Threshold sensitivity | 10/90, 5/95, 20/80; instrumented estimators and native-map summaries | Prepared filtered arrays; no registration or cohort figures |
| Smoothing sensitivity | 4/6/8-mm SPM smoothing; 8/12/16-volume temporal spans; user-data driver | Fixed study acquisition; supplied masks and recorded oxygen/events; regional atlas and tissue/G sensitivity analyses outside scope |
| Convergence/fallback | Per-voxel instrumentation and four-estimator counts | Numerical criterion and fallback encounters, not physiological equilibration |
| Tests/examples | Synthetic arithmetic, isolation, baseline equivalence, SPM and serialization | Do not establish biological accuracy or clinical calibration |
| Study aggregates | LOOCV and sensitivity/convergence aggregates | No participant data or new cohort run |

Read [inputs and outputs](docs/INPUTS_OUTPUTS.md), [statistical definitions](docs/METHODS.md),
[sensitivity and convergence](docs/SENSITIVITY.md), and [validation](docs/VALIDATION.md).
The package contains an internal `SHA256SUMS` for all other files.
Study aggregate tables are provided for [threshold sensitivity](reference/study_threshold_aggregates.csv),
[convergence](reference/study_convergence_aggregates.csv) and [smoothing](reference/study_smoothing_aggregates.csv).
These tables provide cohort-level reference results. Their quantities and
calculation scope are described in the sensitivity documentation.

## Private data

Use [example_private_arrays.m](examples/example_private_arrays.m) only with your
own compatible, aligned data. It never loads bundled study data. For the study,
reproduction requires the study maps from fresh preprocessing, registration/support
provenance, individual tissue masks and native-T1 G masks/scalars. Prepared-array
interfaces do not perform or validate registration. Never mix preprocessing streams.

G-adjusted and conventional deviations answer different model-defined questions.
G is measured from the same subject's response data: LOOCV evaluates conditional
prediction given that measurement. Clinical calibration and biological validity
require independent evaluation.

## Credits and licensing

Code retains its author notices and the [BSD 3-Clause licence](LICENSE).
The external software's citation metadata are retained in [GLIODISSECT_CITATION.cff](GLIODISSECT_CITATION.cff);
they describe GLIODISSECT, not a new authorship claim for this supplement.
Threshold instrumentation and analysis adapters are separately named supplement
implementations. They do not modify the released software or shadow its functions.
