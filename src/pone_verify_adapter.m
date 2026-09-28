function report=pone_verify_adapter()
% Synthetic regression against released and study cores, without private data.
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
old={PONE_DIAGNOSTIC_CONFIG,PONE_DIAGNOSTIC_AUDIT,PONE_DIAGNOSTIC_STATE};c=onCleanup(@()restore(old));
count=0;studyCount=0;
for seed=1:12
 [input,mask]=pone_synthetic_input(seed);
 if seed>6,input.bBOLD(1,1,1,:)=NaN;end
 if seed>9,input.brainLag(2,1,1)=NaN;mask(3,1,1)=false;end
 actual=pone_thresholds(input,mask);
 expected=run_hypoxia_bold_core(input.bBOLD,input.brainLag,input.durationRA,input.initialGlobalLag,input.rowScans,input.shiftedO2,mask);
 fields=fieldnames(expected);
 for k=1:numel(fields)
  assert(isequaln(actual.maps.(fields{k}),expected.(fields{k})),'PONE:Equivalence','Released baseline mismatch: %s',fields{k});count=count+1;
 end
 if seed<=6
  PONE_DIAGNOSTIC_CONFIG=struct('low',.1,'high',.9);PONE_DIAGNOSTIC_AUDIT={};
  study=pone_study_core(input.bBOLD,input.brainLag,input.durationRA,input.initialGlobalLag,input.rowScans,input.shiftedO2,mask);
  fields=fieldnames(actual.maps);
  for k=1:numel(fields)
   if isfield(study,fields{k}),assert(isequaln(actual.maps.(fields{k}),study.(fields{k})),'PONE:StudyEquivalence','Study baseline mismatch: %s',fields{k});studyCount=studyCount+1;end
  end
 end
 assert(all(actual.convergence.finite_output+actual.convergence.NaN_output+actual.convergence.Inf_output==nnz(mask)));
 assert(all(actual.convergence.criterion_satisfied+actual.convergence.criterion_unsatisfied_finite==actual.convergence.evaluable_finite));
 for pair={[.05 .95],[.2 .8]}
  alternative=pone_thresholds(input,mask,pair{1});assert(numel(alternative.audit)==4);
  summary=pone_compare_maps(actual.maps,alternative.maps,mask,'threshold');assert(height(summary)==14);
 end
end
report=struct('synthetic_ensembles',12,'released_field_comparisons',count,'study_core_field_comparisons',studyCount,'exact_equal',true);
end
function restore(old)
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
PONE_DIAGNOSTIC_CONFIG=old{1};PONE_DIAGNOSTIC_AUDIT=old{2};PONE_DIAGNOSTIC_STATE=old{3};
end
