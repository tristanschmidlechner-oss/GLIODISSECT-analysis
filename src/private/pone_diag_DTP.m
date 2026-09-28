% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [DTPmap, stepAmplitude, index10map, index90map, Meanbaseline] = ...
    pone_diag_DTP(BOLD4D, brainlag, DurationRA, lagDiff2, Rowscans)
%PONE_DIAG_DTP Computes 10-90% rise time (DTP) and related maps for each voxel
%
% INPUTS:
%   BOLD4D       - 4D filtered fMRI data (e.g., bBOLD), size: X x Y x Z x T
%   brainlag     - 3D voxel lag map in volumes, same X x Y x Z
%   DurationRA   - Vector [baselineDuration, stepDuration, maxDuration]
%   lagDiff2     - Offset value used in baseline shifting
%   Rowscans     - Time vector of scan TRs (e.g., [1, 2, 3, ..., N])
%
% OUTPUTS:
%   DTPmap           - Map of rise time from 10% to 90% (index90 - index10)
%   stepAmplitude    - BOLD step change amplitude
%   index10map       - Map of 10% response time indices
%   index90map       - Map of 90% response time indices
%   Meanbaseline     - Retained threshold baseline
% This function estimates landmarks for either step. The core uses first-step
% landmarks in refine_lag_step1 to obtain the shared refined total lag.

% Get matrix size
global PONE_DIAGNOSTIC_STATE PONE_DIAGNOSTIC_AUDIT
[MtxX, MtxY, MtxZ, ~] = size(BOLD4D);
auditRows=zeros(MtxX*MtxY*MtxZ,8,'uint16');

% Preallocate outputs
DTPmap = nan(MtxX, MtxY, MtxZ);
stepAmplitude = nan(MtxX, MtxY, MtxZ);
index10map = nan(MtxX, MtxY, MtxZ);
index90map = nan(MtxX, MtxY, MtxZ);
Meanbaseline = nan(MtxX, MtxY, MtxZ);

% Extract durations
Durationbaseline = DurationRA(1);
Durationstepchange = DurationRA(2);
maxDuration = DurationRA(3);

% Loop over all voxels
for ix = 1:MtxX
    for iy = 1:MtxY
        for iz = 1:MtxZ
            av=sub2ind([MtxX MtxY MtxZ],ix,iy,iz);
            PONE_DIAGNOSTIC_STATE=pone_diag_state();
            shift = brainlag(ix, iy, iz);

            % Skip invalid voxels
            if isnan(shift) || shift > maxDuration
                PONE_DIAGNOSTIC_STATE.skipped=true;
                auditRows(av,:)=pone_diag_row(PONE_DIAGNOSTIC_STATE);
                continue
            end

            % Adjust baseline for voxel-specific lag
            Durationbaseline1 = Durationbaseline - lagDiff2 + shift;
            fData = squeeze(BOLD4D(ix, iy, iz, :));

            try
                [finalindex90, finalindex10, finalmeanbaseline, finalstep] = ...
                    pone_diag_crossings(Durationbaseline1, Durationstepchange, Rowscans, fData, shift);

                % Store outputs
                DTPmap(ix, iy, iz) = finalindex90 - finalindex10;
                stepAmplitude(ix, iy, iz) = finalstep;
                index10map(ix, iy, iz) = finalindex10;
                index90map(ix, iy, iz) = finalindex90;
                Meanbaseline(ix, iy, iz) = finalmeanbaseline;
                s=PONE_DIAGNOSTIC_STATE;
                s.criterion=all(isfinite([s.ratioBaseline s.ratioStep])) && all([s.ratioBaseline s.ratioStep]>=0.9999) && all([s.ratioBaseline s.ratioStep]<=1.0001);
                s.bound=s.completed==17 && ~s.criterion;
                s.invalid=~all(isfinite([finalindex10 finalindex90])) || finalindex10<1 || finalindex90<finalindex10 || finalindex90>size(BOLD4D,4) || ~all(isfinite([s.ratioBaseline s.ratioStep]));
                auditRows(av,:)=pone_diag_row(s);

            catch
                PONE_DIAGNOSTIC_STATE.caught=true;
                auditRows(av,:)=pone_diag_row(PONE_DIAGNOSTIC_STATE);
                % Leave NaN on error
                continue
            end
        end
    end
end
PONE_DIAGNOSTIC_AUDIT{end+1}=struct('estimator','DTP','dimensions',[MtxX MtxY MtxZ],'rows',auditRows);
end
