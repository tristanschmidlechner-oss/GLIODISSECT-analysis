function R=pone_thresholds(input,mask,thresholds)
% Prepared filtered data. Supports 10/90, 5/95 and 20/80 thresholds.
if nargin<3,thresholds=[.1 .9];end
assert(size(input.bBOLD,4)==160&&islogical(mask),'PONE:Protocol','160-volume protocol and logical mask required.');
assert(isequal(thresholds,[.1 .9])||isequal(thresholds,[.05 .95])||isequal(thresholds,[.2 .8]),'PONE:Threshold','Use a supported threshold pair.');
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
old={PONE_DIAGNOSTIC_CONFIG,PONE_DIAGNOSTIC_AUDIT,PONE_DIAGNOSTIC_STATE};
cleanup=onCleanup(@()restore(old));
PONE_DIAGNOSTIC_CONFIG=struct('low',thresholds(1),'high',thresholds(2));PONE_DIAGNOSTIC_AUDIT={};
maps=pone_adapter_core(input.bBOLD,input.brainLag,input.durationRA,input.initialGlobalLag,...
 input.rowScans,input.shiftedO2,mask);
R=struct('maps',maps,'audit',{PONE_DIAGNOSTIC_AUDIT},'thresholds',thresholds);
R.convergence=pone_audit_summary(R.audit,maps,mask);
end
function restore(old)
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
PONE_DIAGNOSTIC_CONFIG=old{1};PONE_DIAGNOSTIC_AUDIT=old{2};PONE_DIAGNOSTIC_STATE=old{3};
end
