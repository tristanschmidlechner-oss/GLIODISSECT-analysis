% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [IndexBOLD90b,IndexBOLD10b, meanbaseline,....
    meanstep] = pone_diag_refinement(Startsteptime, IndexBOLD10, IndexBOLD90, Endsteptime, fData,baselineshiftless, shift)
% Bas van Niftrik, Marco Piccirelli, Jorn Fierstra
% Department of Neurosurgery, University Hospital Zurich
%PONE_DIAG_REFINEMENT Update low/high crossings and baseline/step means.
% The baseline ends at the low landmark. The candidate step window runs
% from the high landmark to 45 volumes after the low landmark. Subsequent
% window updates use the fixed 40-volume offset. Crossings use the configured
% low/high fractions and retain the estimator's bounded fallback rules.
% Inputs are one-based scan indices, a voxel time series, and timing offsets.

global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_STATE
PONE_DIAGNOSTIC_STATE.attempts=PONE_DIAGNOSTIC_STATE.attempts+1;
baselinetimeshifted = 1:1:IndexBOLD10;

% create a new array of the stepchange.
steper = IndexBOLD10 + 45; % Candidate step window: 45 volumes after the 10% landmark.

steptimeshift = IndexBOLD90:1:steper;

%% Determine the new mean over the baseline and stepchange. 

Boldbaseline_2 = fData(baselinetimeshifted);
meanboldbaselineold = mean(Boldbaseline_2);

Boldstep_2 = fData(steptimeshift);
meanboldstepold = mean(Boldstep_2);


%% Determine the new shifted 10% and 90%

BOLD10per_2 = (meanboldstepold-meanboldbaselineold).*PONE_DIAGNOSTIC_CONFIG.low;
BOLD10_2c = BOLD10per_2 + meanboldbaselineold;

BOLD90per_2 = (meanboldstepold-meanboldbaselineold).*PONE_DIAGNOSTIC_CONFIG.high;
BOLD90_2c = BOLD90per_2 + meanboldbaselineold;

%% find the new 10% and 90% index value
% Divide between voxels with a positive and a negative signal. 10% can
% either be higher or lower than the baseline.
BOLDh90_2 = -1;
if BOLD90_2c >= meanboldbaselineold
    BOLDh90_2 = find(gt(fData(baselineshiftless:1:steper),BOLD90_2c));
    
elseif BOLD90_2c < meanboldbaselineold
    BOLDh90_2 = find(lt(fData(baselineshiftless:1:steper),BOLD90_2c));
end
%% Determine the new indices for both the 10% and 90%
if ~isempty(BOLDh90_2)
    IndexBOLD90b = BOLDh90_2(1)+ fix(baselineshiftless);
else
    PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(4));
    BOLDh90_2(1) = IndexBOLD90;
    IndexBOLD90b = BOLDh90_2(1)+ fix(baselineshiftless);
end


BOLDh10_2 = -1;
if BOLD10_2c >= meanboldbaselineold
    fData2 = fData(IndexBOLD90b:-1:baselineshiftless);
    BOLDh10_2 =find(lt(fData2,BOLD10_2c));
elseif BOLD10_2c < meanboldbaselineold
    fData2 = fData(IndexBOLD90b:-1:baselineshiftless);
    BOLDh10_2 = find(gt(fData2,BOLD10_2c));
end
if isempty(BOLDh10_2)
    PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(8));
    IndexBOLD10b = fix(baselineshiftless);
else
    IndexBOLD10b = IndexBOLD90b-BOLDh10_2(1);
end


%% create new means
bas = IndexBOLD10b+40;
if (bas-IndexBOLD90b)>= 0
    steptimeshift2 = IndexBOLD90b:1:bas;
else
    PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(16));
    steptimeshift2 = (IndexBOLD90b-5):(IndexBOLD90b+5);
end

fData2 = fData(1:1:IndexBOLD10b);
meanbaseline = mean(fData2);
step_2 = fData(steptimeshift2);
meanstep = mean(step_2);
PONE_DIAGNOSTIC_STATE.completed=PONE_DIAGNOSTIC_STATE.completed+1;
PONE_DIAGNOSTIC_STATE.ratioBaseline=meanbaseline/PONE_DIAGNOSTIC_STATE.lastBaseline;
PONE_DIAGNOSTIC_STATE.ratioStep=meanstep/PONE_DIAGNOSTIC_STATE.lastStep;
PONE_DIAGNOSTIC_STATE.lastBaseline=meanbaseline;PONE_DIAGNOSTIC_STATE.lastStep=meanstep;
end
