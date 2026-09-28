% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function G = pone_g_scalar(steadyStateO2NormT1, greyWhiteMaskT1)
%PONE_G_SCALAR Median global O2-normalised steady-state response.

validateattributes(greyWhiteMaskT1, {'logical','numeric'}, {'size', size(steadyStateO2NormT1)});
values = steadyStateO2NormT1(logical(greyWhiteMaskT1) & isfinite(steadyStateO2NormT1));
if isempty(values)
    error('PONE:EmptyGMask', 'No finite values exist inside the G-scalar mask.');
end
G = median(double(values));
end
