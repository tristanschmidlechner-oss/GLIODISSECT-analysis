% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function maps = pone_study_core(bBOLD, brainLag, durationRA, ...
    initialGlobalLag, rowScans, shiftedO2, analysisMask)
%PONE_STUDY_CORE Compute temporal and amplitude maps using study conventions.
%
% This function is cohort-independent. Healthy-control and patient runners
% construct their masks separately and then call this shared numerical core.
%
% Inputs
%   bBOLD            Filtered four-dimensional BOLD data (X-by-Y-by-Z-by-T)
%   brainLag         Voxelwise initial lag map in volumes, before the
%                    +0.001/-0.001 indexing arithmetic (GLIODISSECT LAG_DEFINITIONS)
%   durationRA       Five protocol phase durations in volumes
%   initialGlobalLag Initial global lag used for the first DTP estimate
%   rowScans         Scan indices, normally 1:T
%   shiftedO2        End-tidal O2 sampled at the BOLD acquisition times
%   analysisMask     Logical EPI-space cohort-specific analysis mask
%
% Output
%   maps             Structure containing DTP, DTR, first-step lag,
%                    steady-state and whole-stimulus-average maps.

validateattributes(bBOLD, {'numeric'}, {'nonempty', 'real'}, mfilename, 'bBOLD');
spatialSize = [size(bBOLD, 1), size(bBOLD, 2), size(bBOLD, 3)];
validateattributes(brainLag, {'numeric'}, {'size', spatialSize}, mfilename, 'brainLag');
validateattributes(durationRA, {'numeric'}, {'vector', 'numel', 5, 'finite'}, mfilename, 'durationRA');
validateattributes(rowScans, {'numeric'}, {'vector', 'numel', size(bBOLD, 4)}, mfilename, 'rowScans');
if ~isa(shiftedO2, 'function_handle')
    validateattributes(shiftedO2, {'numeric'}, {'vector', 'numel', size(bBOLD, 4)}, mfilename, 'shiftedO2');
end
validateattributes(analysisMask, {'logical', 'numeric'}, {'size', size(brainLag)}, mfilename, 'analysisMask');

analysisMask = logical(analysisMask);
expandedMask = repmat(analysisMask, 1, 1, 1, size(bBOLD, 4));
bBOLD(~expandedMask) = NaN;
initialLag = brainLag;

% Ordered floating-point offset arithmetic is part of the indexing contract.
% It can affect integer boundaries; see GLIODISSECT docs/LAG_DEFINITIONS.md.
brainLag = brainLag + 0.001;
brainLag = brainLag - 0.001;

[dtp1, stepMean1, index10_1, index90_1, baselineMean1] = ...
    pone_diag_DTP(bBOLD, brainLag, durationRA, initialGlobalLag, rowScans);
refinedLagRaw = refine_lag_step1(brainLag, index10_1, durationRA(1), initialGlobalLag);
% Calculate the nonnegative QC adjustment before applying offset arithmetic.
clippedLagAdjustment = nan(size(brainLag));
validRefinement = isfinite(refinedLagRaw);
clippedLagAdjustment(validRefinement) = max(brainLag(validRefinement) - refinedLagRaw(validRefinement), 0);
refinedLag = refinedLagRaw + 0.001;
refinedLag = refinedLag - 0.001;
clippedLagAdjustment = clippedLagAdjustment + 0.001;
clippedLagAdjustment = clippedLagAdjustment - 0.001;
signedLagAdjustment = brainLag - refinedLag;
validRefinedLag = analysisMask & isfinite(refinedLag);
if ~any(validRefinedLag(:))
    error('PONE:NoValidDelay', 'No finite first-step refined lags exist inside the analysis mask.');
end
globalDelay = round(prctile(refinedLag(validRefinedLag), 5));
globalDelay = max(0, min(globalDelay, size(bBOLD, 4) - 5));
validArrival = analysisMask & isfinite(index10_1) & isfinite(refinedLag);
referenceI10 = prctile(index10_1(validArrival), 5);
relativeArrivalDelay = index10_1 - referenceI10;
relativeArrivalDelay(~validArrival) = NaN;
% Resample the recorded O2 trace after the refined global lag is available.
if isa(shiftedO2, 'function_handle')
    shiftedO2 = shiftedO2(globalDelay);
    validateattributes(shiftedO2, {'numeric'}, {'vector', 'numel', size(bBOLD, 4)});
end

adjustedDuration = durationRA(:)';
adjustedDuration(1) = adjustedDuration(1) + round(globalDelay - initialGlobalLag);
durationStep1 = adjustedDuration(1:3);
durationStep2 = [sum(adjustedDuration(1:3)), adjustedDuration(4), adjustedDuration(5)];

[dtr1, return90_1, return10_1, postBaselineMean1] = pone_diag_DTR( ...
    bBOLD, refinedLag, durationStep1, globalDelay, rowScans, index90_1, stepMean1);
[dtp2, stepMean2, index10_2, index90_2, baselineMean2] = ...
    pone_diag_DTP(bBOLD, refinedLag, durationStep2, globalDelay, rowScans);
[dtr2, return90_2, return10_2, postBaselineMean2] = pone_diag_DTR( ...
    bBOLD, refinedLag, durationStep2, globalDelay, rowScans, index90_2, stepMean2);

[steadyState1, steadyState2, steadyStateDuration1, steadyStateDuration2, ...
    steadyStateO2Norm1, steadyStateO2Norm2, deltaPetO2_1, deltaPetO2_2] = ...
    OxygenDoubleStep_SS(bBOLD, index90_1, return90_1, baselineMean1, ...
    index90_2, return90_2, refinedLag, durationStep1, durationStep2, ...
    globalDelay, shiftedO2);

[stimulusAverage1, stimulusAverage2] = compute_stimulus_average( ...
    refinedLag, adjustedDuration, globalDelay, rowScans, bBOLD);

maps = struct();
maps.DTP_step1 = apply_mask(dtp1, analysisMask);
% Study-core DTP2 is returned without an output mask; see docs/VALIDATION.md.
maps.DTP_step2 = dtp2;
maps.DTR_step1 = apply_mask(dtr1, analysisMask);
maps.DTR_step2 = apply_mask(dtr2, analysisMask);
maps.initial_lag = apply_mask(initialLag, analysisMask);
maps.refined_lag_step1 = apply_mask(refinedLag, analysisMask);
maps.signed_lag_adjustment_step1 = apply_mask(signedLagAdjustment, analysisMask);
maps.clipped_lag_adjustment_step1 = apply_mask(clippedLagAdjustment, analysisMask);
maps.relative_arrival_delay_step1 = relativeArrivalDelay;
maps.relative_arrival_reference_i10_step1 = referenceI10;
maps.steady_state_step1 = apply_mask(steadyState1, analysisMask);
maps.steady_state_step2 = apply_mask(steadyState2, analysisMask);
maps.steady_state_O2norm_step1 = apply_mask(steadyStateO2Norm1, analysisMask);
maps.steady_state_O2norm_step2 = apply_mask(steadyStateO2Norm2, analysisMask);
maps.steady_state_duration_step1 = apply_mask(steadyStateDuration1, analysisMask);
maps.steady_state_duration_step2 = apply_mask(steadyStateDuration2, analysisMask);
maps.stimulus_average_step1 = apply_mask(stimulusAverage1, analysisMask);
maps.stimulus_average_step2 = apply_mask(stimulusAverage2, analysisMask);
maps.stimulus_average_O2norm_step1 = apply_mask(stimulusAverage1 / deltaPetO2_1, analysisMask);
maps.stimulus_average_O2norm_step2 = apply_mask(stimulusAverage2 / deltaPetO2_2, analysisMask);
maps.index10_DTP_step1 = apply_mask(index10_1, analysisMask);
maps.index90_DTP_step1 = apply_mask(index90_1, analysisMask);
maps.index10_DTP_step2 = apply_mask(index10_2, analysisMask);
maps.index90_DTP_step2 = apply_mask(index90_2, analysisMask);
maps.index90_DTR_step1 = apply_mask(return90_1, analysisMask);
maps.index10_DTR_step1 = apply_mask(return10_1, analysisMask);
maps.index90_DTR_step2 = apply_mask(return90_2, analysisMask);
maps.index10_DTR_step2 = apply_mask(return10_2, analysisMask);
maps.baseline_mean_step1 = apply_mask(baselineMean1, analysisMask);
maps.baseline_mean_step2 = apply_mask(baselineMean2, analysisMask);
maps.post_baseline_mean_step1 = apply_mask(postBaselineMean1, analysisMask);
maps.post_baseline_mean_step2 = apply_mask(postBaselineMean2, analysisMask);
maps.delta_PetO2_step1 = deltaPetO2_1;
maps.delta_PetO2_step2 = deltaPetO2_2;
maps.global_refined_lag_step1_volumes = globalDelay;
maps.analysis_mask = analysisMask;
end

function image = apply_mask(image, mask)
image(~mask) = NaN;
end
