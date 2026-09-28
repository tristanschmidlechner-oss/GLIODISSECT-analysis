# Statistical definitions

At each eligible voxel OLS fits Y = beta0 + beta1*G from valid observations.
The double-precision sufficient statistics use centered sums cxx and cxy,
beta1=cxy/cxx, beta0=mean(Y)-beta1*mean(G), SSE=max(SST-beta1*cxy,0),
residual SD=sqrt(SSE/(n-2)). Positive cxx and finite coefficients are required.
The fixed floor is MATLAB prctile(positive finite eligible residual SD,5).
Constant response yields undefined R²; zero/negative Q² is not suppressed.
The conventional model uses single-precision mean and sample std(Y,0,2), with
NaNs omitted and a P5 floor from its own eligible positive sample SDs.

For the study's 20 LOOCV folds, select the 19 training controls before computing
anatomical coverage (>=16), valid response coverage (>=16), OLS or residual floor.
Held-out G is used only for prediction; held-out response/anatomy only for scoring.
The comparator is the mean of the exact same valid training observations used
by OLS at that voxel. Prediction support is training eligibility intersected with
the held-out individual tissue and finite-response support. The full-control fit
and its floor never enter a fold.

Subject RMSE = sqrt(mean(error²)) on that fold's scored support. Pooled RMSE =
sqrt(sum of all scored squared errors / total scored predictions), including
predictions at voxels that do not meet Q² coverage. It is not mean subject RMSE.
Reduction = 100*(RMSE_mean-RMSE_G)/RMSE_mean. Improvement requires strictly lower RMSE.
Voxel Q²=1-PRESS/TSS about the mean of eligible held-out observations, requiring
>=17 predictions and positive TSS. The mean comparator uses the same denominator.
The residual floor affects standardized errors and descriptive |z|>2 rates, not
predictions, Q² or RMSE. Missing predictions are not filled with zeros.

In-sample R² is computed separately using all 20 controls, with anatomy and
response coverage >=16/20. It is not cross-validated Q². Other cohort sizes retain
these definitions with explicitly supplied integer counts; this is a software
generalization, not cohort validation of those choices.

The [saved study aggregate reference](../reference/study_loocv_aggregates.json)
records the study results: Q² median 0.345005146746427; RMSE reduction
16.0264138360177%; 19/20 controls improved; in-sample R² median
0.490301917600222; pooled RMSE G 0.0184303259354605 and mean
0.0219477656932145. These results describe conditional prediction given each subject's G scalar.
