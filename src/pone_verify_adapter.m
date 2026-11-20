function report=pone_verify_adapter()
% Synthetic regression against GLIODISSECT v1.0.0, without private data.
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
old={PONE_DIAGNOSTIC_CONFIG,PONE_DIAGNOSTIC_AUDIT,PONE_DIAGNOSTIC_STATE};c=onCleanup(@()restore(old));
count=0;averageCount=0;oxygenCount=0;
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
 averageCount=averageCount+verify_averages(input,actual.maps,expected);
 assert(all(actual.convergence.finite_output+actual.convergence.NaN_output+actual.convergence.Inf_output==nnz(mask)));
 assert(all(actual.convergence.criterion_satisfied+actual.convergence.criterion_unsatisfied_finite==actual.convergence.evaluable_finite));
 for pair={[.05 .95],[.2 .8]}
  alternative=pone_thresholds(input,mask,pair{1});assert(numel(alternative.audit)==4);
  summary=pone_compare_maps(actual.maps,alternative.maps,mask,'threshold');assert(height(summary)==14);
 end
end
% Exercise zero and reversed oxygen differences, and the callback interface.
for scenario=1:3
 [input,mask]=pone_synthetic_input(5);
 oxygen=input.shiftedO2;
 if scenario==1,input.shiftedO2=100*ones(size(oxygen));end
 if scenario==2,input.shiftedO2=200-oxygen;end
 if scenario==3,input.shiftedO2=@(~)oxygen;end
 actual=pone_thresholds(input,mask);
 expected=run_hypoxia_bold_core(input.bBOLD,input.brainLag,input.durationRA,input.initialGlobalLag,input.rowScans,input.shiftedO2,mask);
 fields=fieldnames(expected);
 for k=1:numel(fields)
  assert(isequaln(actual.maps.(fields{k}),expected.(fields{k})),'PONE:Equivalence','Released oxygen-case mismatch: %s',fields{k});oxygenCount=oxygenCount+1;
 end
 if scenario==1
  assert(expected.delta_PetO2_step1==0&&expected.delta_PetO2_step2==0);
  assert(all(isnan(actual.maps.stimulus_average_O2norm_step1(:)))&&all(isnan(actual.maps.stimulus_average_O2norm_step2(:))));
 elseif scenario==2
  assert(expected.delta_PetO2_step1<0&&expected.delta_PetO2_step2<0);
 end
 averageCount=averageCount+verify_averages(input,actual.maps,expected);
end
report=struct('synthetic_ensembles',12,'released_field_comparisons',count,...
 'oxygen_cases',3,'oxygen_released_field_comparisons',oxygenCount,...
 'whole_average_field_checks',averageCount,'whole_average_absolute_tolerance',1e-11,'exact_equal',true);
end

function count=verify_averages(input,actual,released)
% Independent phase-weight oracle. Timing/support and oxygen differences come
% from the released implementation; the supplement's averaging helper is not used.
% These fixtures have complete, nonempty phases inside the acquisition.
raw=nan([size(released.validity_mask),2]);
duration=input.durationRA;
globalLag=released.global_refined_lag_step1_volumes;
duration(1)=duration(1)+round(globalLag-input.initialGlobalLag);
for v=find(released.validity_mask(:))'
 [x,y,z]=ind2sub(size(released.validity_mask),v);
 first=duration(1)-globalLag+released.refined_lag_step1(v);
 edges=[0 floor(first+cumsum([0 duration(2:4)])) size(input.bBOLD,4)];
 assert(all(diff(edges)>0),'PONE:OracleFixture','Oracle requires five complete phases.');
 membership=(input.rowScans(:)>edges(1:5))&(input.rowScans(:)<=edges(2:6));
 weights=double(membership)./sum(membership,1);
 phaseMeans=double(reshape(input.bBOLD(x,y,z,:),[],1))'*weights;
 baseline=phaseMeans([1 3 5])*ones(3,1)/3;
 assert(isfinite(baseline)&&baseline~=0,'PONE:OracleFixture','Finite nonzero baseline required.');
 raw(x,y,z,1)=100*(phaseMeans(2)-baseline)/baseline;
 raw(x,y,z,2)=100*(phaseMeans(4)-baseline)/baseline;
end
assert(any(isfinite(raw(:))),'PONE:OracleFixture','Fixture must exercise finite averages.');
names={'stimulus_average_step1','stimulus_average_step2',...
 'stimulus_average_O2norm_step1','stimulus_average_O2norm_step2'};
values={raw(:,:,:,1),raw(:,:,:,2),raw(:,:,:,1)/released.delta_PetO2_step1,raw(:,:,:,2)/released.delta_PetO2_step2};
for k=1:numel(names)
 expected=values{k};expected(~released.validity_mask|~isfinite(expected))=NaN;
 observed=actual.(names{k});
 assert(isequal(size(observed),size(expected))&&isequal(isnan(observed),isnan(expected)),...
  'PONE:AverageSupport','Whole-average support mismatch: %s',names{k});
 finite=isfinite(expected);
 assert(all(abs(observed(finite)-expected(finite))<1e-11),...
  'PONE:AverageOracle','Whole-average arithmetic mismatch: %s',names{k});
end
count=numel(names);
end
function restore(old)
global PONE_DIAGNOSTIC_CONFIG PONE_DIAGNOSTIC_AUDIT PONE_DIAGNOSTIC_STATE
PONE_DIAGNOSTIC_CONFIG=old{1};PONE_DIAGNOSTIC_AUDIT=old{2};PONE_DIAGNOSTIC_STATE=old{3};
end
