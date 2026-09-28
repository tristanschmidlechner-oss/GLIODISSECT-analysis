% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [stepAverage1, stepAverage2] = compute_stimulus_average( ...
    voxelDelay, durationRA, globalDelay, rowScans, bBOLD)
%COMPUTE_STIMULUS_AVERAGE Whole-stimulus response averages.
%
% Baseline signal is the equally weighted mean of the three baseline-period
% means; each response is the percentage change during its full step.

spatialSize = [size(bBOLD,1), size(bBOLD,2), size(bBOLD,3)];
stepAverage1 = nan(spatialSize);
stepAverage2 = nan(spatialSize);
numberOfScans = size(bBOLD,4);

for ix = 1:spatialSize(1)
    for iy = 1:spatialSize(2)
        for iz = 1:spatialSize(3)
            delay = voxelDelay(ix,iy,iz);
            if ~isfinite(delay); continue; end
            baseline1Duration = durationRA(1) - globalDelay + delay;
            step1Duration = durationRA(2);
            baseline2Duration = durationRA(3);
            step2Duration = durationRA(4);

            endpoints = [baseline1Duration - 1, ...
                baseline1Duration, baseline1Duration + step1Duration - 1, ...
                baseline1Duration + step1Duration, ...
                baseline1Duration + step1Duration + baseline2Duration - 1, ...
                baseline1Duration + step1Duration + baseline2Duration, ...
                baseline1Duration + step1Duration + baseline2Duration + step2Duration - 1, ...
                baseline1Duration + step1Duration + baseline2Duration + step2Duration];
            index = arrayfun(@(value) first_greater(rowScans, value, numberOfScans), endpoints);
            ranges = {1:index(1), index(2):index(3), index(4):index(5), ...
                index(6):index(7), index(8):numberOfScans};
            if any(cellfun(@isempty, ranges)); continue; end

            series = squeeze(bBOLD(ix,iy,iz,:));
            baselineMean = mean([mean(series(ranges{1})), ...
                mean(series(ranges{3})), mean(series(ranges{5}))]);
            if ~isfinite(baselineMean) || baselineMean == 0; continue; end
            stepAverage1(ix,iy,iz) = (mean(series(ranges{2})) - baselineMean) / baselineMean * 100;
            stepAverage2(ix,iy,iz) = (mean(series(ranges{4})) - baselineMean) / baselineMean * 100;
        end
    end
end
end

function index = first_greater(rowScans, value, numberOfScans)
index = find(rowScans > value, 1, 'first');
if isempty(index); index = numberOfScans; end
index = min(max(index,1),numberOfScans);
end
