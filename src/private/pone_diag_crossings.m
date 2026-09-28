% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [finalindex90, finalindex10,...
    finalmeanbaseline, finalmeanstep] = pone_diag_crossings(Durationbaseline1,...
   Durationstepchange, Rowscans, fData, shift) 
% Bas van Niftrik, Marco Piccirelli, Jorn Fierstra
% Department of Neurosurgery, University Hospital Zurich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate threshold crossings with bounded baseline and step refinement.
% Window boundaries use the first matching sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter input
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_STATE
Durationbaseline  = Durationbaseline1;
DurationStepchange = Durationstepchange;

%% calculate durations of each step
Startbaseline1time = find(gt(Rowscans, shift));
Startbaseline1time = Startbaseline1time(1);

Endbaseline1time   = find(gt(Rowscans, (Durationbaseline-1)));
Endbaseline1time   = Endbaseline1time(1);
baseline1timeindex = Startbaseline1time:Endbaseline1time;

Startsteptime = find(gt(Rowscans, Durationbaseline));
Startsteptime = Startsteptime(1);          

Endsteptime   = find(gt(Rowscans, (Durationbaseline+DurationStepchange+5)));
Endsteptime   = Endsteptime(1);
Steptime      = Startsteptime:Endsteptime;

baselineshiftless = fix(Durationbaseline)-fix(shift);

%% determine the BOLD during the first two phases
BOLDbaseline1    = fData(baseline1timeindex);
MeanBOLDbaseline1 = mean(BOLDbaseline1);

BOLDstep     = fData(Steptime);
MeanBOLDStep = mean(BOLDstep);

%% Calculate the 10% and 90% value
BOLD10per = (MeanBOLDStep-MeanBOLDbaseline1).*PONE_DIAGNOSTIC_CONFIG.low;
BOLD10    = BOLD10per + MeanBOLDbaseline1;

BOLD90per = (MeanBOLDStep-MeanBOLDbaseline1).*PONE_DIAGNOSTIC_CONFIG.high;
BOLD90    = BOLD90per + MeanBOLDbaseline1;

%% find the 90% index value
BOLDh90 = -1.0;
if BOLD90 >= MeanBOLDbaseline1
    BOLDh90 = find(gt(fData(baselineshiftless:Endsteptime),BOLD90));
elseif BOLD90 < MeanBOLDbaseline1
    BOLDh90 = find(lt(fData(baselineshiftless:Endsteptime),BOLD90));
end
BOLDh90_10 = BOLDh90;

if ~isempty(BOLDh90_10)
    IndexBOLD90 = BOLDh90_10(1) + baselineshiftless;
else
    PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(1));
        BOLDh90_10(1) = size(Steptime,2)-5;
    IndexBOLD90   = BOLDh90_10(1) + baselineshiftless;
end

BOLDh10    = -1;
IndexBOLD10 = 1;

%% find the 10% index value
if BOLD10 >= MeanBOLDbaseline1
    BOLDh10 = find(lt(fData(IndexBOLD90:-1:baselineshiftless),BOLD10));
    if isempty(BOLDh10)
        PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(2));
        IndexBOLD10 = fix(baselineshiftless);
    else
        IndexBOLD10 = IndexBOLD90 - BOLDh10(1) + 1;
    end
elseif BOLD10 < MeanBOLDbaseline1
    BOLDh10 = find(gt(fData(IndexBOLD90:-1:baselineshiftless),BOLD10));
    if isempty(BOLDh10)
        PONE_DIAGNOSTIC_STATE.fallback=bitor(PONE_DIAGNOSTIC_STATE.fallback,uint16(2));
        IndexBOLD10 = fix(baselineshiftless);
    else
        IndexBOLD10 = IndexBOLD90 - BOLDh10(1) - 1;
    end
end

%% Two fixed refinements followed by up to 15 additional refinements
[IndexBOLD90b,IndexBOLD10b,meanboldbaseline2,meanboldstep2] = ...
    pone_diag_refinement(Startsteptime,IndexBOLD10,IndexBOLD90,Endsteptime,fData,baselineshiftless,shift);

[IndexBOLD90c,IndexBOLD10c,meanboldbaseline3,meanboldstep3] = ...
    pone_diag_refinement(Startsteptime,IndexBOLD10b,IndexBOLD90b,Endsteptime,fData,baselineshiftless,shift);

baselineDiff = meanboldbaseline3./meanboldbaseline2;
if (baselineDiff>1.0001||baselineDiff<0.9999) || ((meanboldstep3./meanboldstep2)>1.0001||(meanboldstep3./meanboldstep2)<0.9999)
    [IndexBOLD90d,IndexBOLD10d,meanboldbaseline4,meanboldstep4] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10c,IndexBOLD90c,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90d=IndexBOLD90c; meanboldstep4=meanboldstep3;
    IndexBOLD10d=IndexBOLD10c; meanboldbaseline4=meanboldbaseline3;
end

baselineDiff2 = meanboldbaseline4./meanboldbaseline3;
if (baselineDiff2>1.0001||baselineDiff2<0.9999) || ((meanboldstep4./meanboldstep3)>1.0001||(meanboldstep4./meanboldstep3)<0.9999)
    [IndexBOLD90e,IndexBOLD10e,meanboldbaseline5,meanboldstep5] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10d,IndexBOLD90d,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90e=IndexBOLD90d; meanboldstep5=meanboldstep4;
    IndexBOLD10e=IndexBOLD10d; meanboldbaseline5=meanboldbaseline4;
end

baselineDiff3 = meanboldbaseline5./meanboldbaseline4;
if (baselineDiff3>1.0001||baselineDiff3<0.9999) || ((meanboldstep5./meanboldstep4)>1.0001||(meanboldstep5./meanboldstep4)<0.9999)
    [IndexBOLD90f,IndexBOLD10f,meanboldbaseline6,meanboldstep6] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10e,IndexBOLD90e,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90f=IndexBOLD90e; meanboldstep6=meanboldstep5;
    IndexBOLD10f=IndexBOLD10e; meanboldbaseline6=meanboldbaseline5;
end

baselineDiff4 = meanboldbaseline6./meanboldbaseline5;
if (baselineDiff4>1.0001||baselineDiff4<0.9999) || ((meanboldstep6./meanboldstep5)>1.0001||(meanboldstep6./meanboldstep5)<0.9999)
    [IndexBOLD90g,IndexBOLD10g,meanboldbaseline7,meanboldstep7] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10f,IndexBOLD90f,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90g=IndexBOLD90f; meanboldstep7=meanboldstep6;
    IndexBOLD10g=IndexBOLD10f; meanboldbaseline7=meanboldbaseline6;
end

baselineDiff5 = meanboldbaseline7./meanboldbaseline6;
if (baselineDiff5>1.0001||baselineDiff5<0.9999) || ((meanboldstep7./meanboldstep6)>1.0001||(meanboldstep7./meanboldstep6)<0.9999)
    [IndexBOLD90h,IndexBOLD10h,meanboldbaseline8,meanboldstep8] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10g,IndexBOLD90g,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90h=IndexBOLD90g; meanboldstep8=meanboldstep7;
    IndexBOLD10h=IndexBOLD10g; meanboldbaseline8=meanboldbaseline7;
end

baselineDiff6 = meanboldbaseline8./meanboldbaseline7;
if (baselineDiff6>1.0001||baselineDiff6<0.9999) || ((meanboldstep8./meanboldstep7)>1.0001||(meanboldstep8./meanboldstep7)<0.9999)
    [IndexBOLD90i,IndexBOLD10i,meanboldbaseline9,meanboldstep9] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10h,IndexBOLD90h,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90i=IndexBOLD90h; meanboldstep9=meanboldstep8;
    IndexBOLD10i=IndexBOLD10h; meanboldbaseline9=meanboldbaseline8;
end

baselineDiff7 = meanboldbaseline9./meanboldbaseline8;
if (baselineDiff7>1.0001||baselineDiff7<0.9999) || ((meanboldstep9./meanboldstep8)>1.0001||(meanboldstep9./meanboldstep8)<0.9999)
    [IndexBOLD90j,IndexBOLD10j,meanboldbaseline10,meanboldstep10] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10i,IndexBOLD90i,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90j=IndexBOLD90i; meanboldstep10=meanboldstep9;
    IndexBOLD10j=IndexBOLD10i; meanboldbaseline10=meanboldbaseline9;
end

baselineDiff8 = meanboldbaseline10./meanboldbaseline9;
if (baselineDiff8>1.0001||baselineDiff8<0.9999) || ((meanboldstep10./meanboldstep9)>1.0001||(meanboldstep10./meanboldstep9)<0.9999)
    [IndexBOLD90k,IndexBOLD10k,meanboldbaseline11,meanboldstep11] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10j,IndexBOLD90j,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90k=IndexBOLD90j; meanboldstep11=meanboldstep10;
    IndexBOLD10k=IndexBOLD10j; meanboldbaseline11=meanboldbaseline10;
end

baselineDiff9 = meanboldbaseline11./meanboldbaseline10;
if (baselineDiff9>1.0001||baselineDiff9<0.9999) || ((meanboldstep11./meanboldstep10)>1.0001||(meanboldstep11./meanboldstep10)<0.9999)
    [IndexBOLD90l,IndexBOLD10l,meanboldbaseline12,meanboldstep12] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10k,IndexBOLD90k,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90l=IndexBOLD90k; meanboldstep12=meanboldstep11;
    IndexBOLD10l=IndexBOLD10k; meanboldbaseline12=meanboldbaseline11;
end

baselineDiff10 = meanboldbaseline12./meanboldbaseline11;
if (baselineDiff10>1.0001||baselineDiff10<0.9999) || ((meanboldstep12./meanboldstep11)>1.0001||(meanboldstep12./meanboldstep11)<0.9999)
    [IndexBOLD90m,IndexBOLD10m,meanboldbaseline13,meanboldstep13] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10l,IndexBOLD90l,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90m=IndexBOLD90l; meanboldstep13=meanboldstep12;
    IndexBOLD10m=IndexBOLD10l; meanboldbaseline13=meanboldbaseline12;
end

baselineDiff11 = meanboldbaseline13./meanboldbaseline12;
if (baselineDiff11>1.0001||baselineDiff11<0.9999) || ((meanboldstep13./meanboldstep12)>1.0001||(meanboldstep13./meanboldstep12)<0.9999)
    [IndexBOLD90n,IndexBOLD10n,meanboldbaseline14,meanboldstep14] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10m,IndexBOLD90m,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90n=IndexBOLD90m; meanboldstep14=meanboldstep13;
    IndexBOLD10n=IndexBOLD10m; meanboldbaseline14=meanboldbaseline13;
end

baselineDiff12 = meanboldbaseline14./meanboldbaseline13;
if (baselineDiff12>1.0001||baselineDiff12<0.9999) || ((meanboldstep14./meanboldstep13)>1.0001||(meanboldstep14./meanboldstep13)<0.9999)
    [IndexBOLD90o,IndexBOLD10o,meanboldbaseline15,meanboldstep15] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10n,IndexBOLD90n,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90o=IndexBOLD90n; meanboldstep15=meanboldstep14;
    IndexBOLD10o=IndexBOLD10n; meanboldbaseline15=meanboldbaseline14;
end

baselineDiff13 = meanboldbaseline15./meanboldbaseline14;
if (baselineDiff13>1.0001||baselineDiff13<0.9999) || ((meanboldstep15./meanboldstep14)>1.0001||(meanboldstep15./meanboldstep14)<0.9999)
    [IndexBOLD90p,IndexBOLD10p,meanboldbaseline16,meanboldstep16] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10o,IndexBOLD90o,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90p=IndexBOLD90o; meanboldstep16=meanboldstep15;
    IndexBOLD10p=IndexBOLD10o; meanboldbaseline16=meanboldbaseline15;
end

baselineDiff14 = meanboldbaseline16./meanboldbaseline15;
if (baselineDiff14>1.0001||baselineDiff14<0.9999) || ((meanboldstep16./meanboldstep15)>1.0001||(meanboldstep16./meanboldstep15)<0.9999)
    [IndexBOLD90q,IndexBOLD10q,meanboldbaseline17,meanboldstep17] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10p,IndexBOLD90p,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90q=IndexBOLD90p; meanboldstep17=meanboldstep16;
    IndexBOLD10q=IndexBOLD10p; meanboldbaseline17=meanboldbaseline16;
end

baselineDiff15 = meanboldbaseline17./meanboldbaseline16;
if (baselineDiff15>1.0001||baselineDiff15<0.9999) || ((meanboldstep17./meanboldstep16)>1.0001||(meanboldstep17./meanboldstep16)<0.9999)
    [IndexBOLD90r,IndexBOLD10r,meanboldbaseline18,meanboldstep18] = ...
        pone_diag_refinement(Startsteptime,IndexBOLD10q,IndexBOLD90q,Endsteptime,fData,baselineshiftless,shift);
else
    IndexBOLD90r=IndexBOLD90q; meanboldstep18=meanboldstep17;
    IndexBOLD10r=IndexBOLD10q; meanboldbaseline18=meanboldbaseline17;
end

%% assign outputs
finalmeanbaseline = meanboldbaseline18;
finalmeanstep     = meanboldstep18;
finalindex90      = IndexBOLD90r;
finalindex10      = IndexBOLD10r;
end
