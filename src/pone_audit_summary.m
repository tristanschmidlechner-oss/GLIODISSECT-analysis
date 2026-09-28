function T=pone_audit_summary(audit,maps,mask)
% Convergence/fallback counters; finite output is not physiological convergence.
assert(numel(audit)==4,'PONE:Audit','Four estimator records required.');
estimator={'DTP1','DTR1','DTP2','DTR2'};mapfield={'DTP_step1','DTR_step1','DTP_step2','DTR_step2'};
A=struct([]);
for e=1:4
  q=double(audit{e}.rows);assert(size(q,1)==numel(mask));vals=maps.(mapfield{e});dom=mask(:);attempt=q(:,1)>0&dom;finite=isfinite(vals(:))&dom;evaluable=finite&~logical(q(:,7))&~logical(q(:,8));criterion=logical(q(:,3))&evaluable;bound=17;if contains(estimator{e},'DTR'),bound=15;end
  unsat=evaluable&~criterion;reached=q(:,2)==bound&dom;d=q(attempt,2);pr=prctile(d,[0 25 50 75 100]);
  r=struct('estimator',string(estimator{e}),'mask_voxels',nnz(dom),'attempted',nnz(attempt),'evaluable_finite',nnz(evaluable),'refinement_min',pr(1),'refinement_p25',pr(2),'refinement_median',pr(3),'refinement_p75',pr(4),'refinement_IQR',pr(4)-pr(2),'refinement_max',pr(5),'criterion_satisfied',nnz(criterion),'criterion_percent_evaluable',100*nnz(criterion)/nnz(evaluable),'criterion_unsatisfied_finite',nnz(unsat),'criterion_unsatisfied_percent_evaluable',100*nnz(unsat)/nnz(evaluable),'reached_bound',nnz(reached),'reached_bound_percent_mask',100*nnz(reached)/nnz(dom),'bounded_unsatisfied',nnz(unsat&reached),'any_fallback',nnz(q(:,5)>0&dom),'fallback_percent_mask',100*nnz(q(:,5)>0&dom)/nnz(dom),'invalid_crossing_or_nonfinite_criterion',nnz(q(:,6)>0&dom),'caught_error',nnz(q(:,7)>0&dom),'skipped_invalid_input_or_window',nnz(q(:,8)>0&dom),'finite_output',nnz(finite),'zero_output',nnz(vals(:)==0&dom),'NaN_output',nnz(isnan(vals(:))&dom),'Inf_output',nnz(isinf(vals(:))&dom));
  if e==1,a10=maps.index10_DTP_step1;a90=maps.index90_DTP_step1;elseif e==2,a10=maps.index10_DTR_step1;a90=maps.index90_DTR_step1;elseif e==3,a10=maps.index10_DTP_step2;a90=maps.index90_DTP_step2;else,a10=maps.index10_DTR_step2;a90=maps.index90_DTR_step2;end;bad=~isfinite(a10(:))|~isfinite(a90(:))|a10(:)<1|a90(:)<1|a10(:)>160|a90(:)>160;if contains(estimator{e},'DTP'),bad=bad|a90(:)<a10(:);else,bad=bad|a10(:)<a90(:);end;r.missing_or_invalid_final_crossing=nnz(dom&bad);r.internal_invalid_flag_with_valid_final_indices=nnz(dom&~bad&q(:,6)>0);
  for b=0:4,r.(sprintf('fallback_bit%d_voxels',b))=nnz(bitget(uint16(q(:,5)),b+1)&dom);end
  crossbits=uint16(15);windowbits=uint16(16);if contains(estimator{e},'DTR'),crossbits=uint16(19);windowbits=uint16(12);end;r.crossing_fallback=nnz(dom&bitand(uint16(q(:,5)),crossbits)>0);r.window_rule_fallback=nnz(dom&bitand(uint16(q(:,5)),windowbits)>0);r.crossing_fallback_percent_mask=100*r.crossing_fallback/nnz(dom);assert(r.finite_output+r.NaN_output+r.Inf_output==nnz(dom));

 if isempty(A),A=r;else,A(end+1)=r;end
end
T=struct2table(A);
end
