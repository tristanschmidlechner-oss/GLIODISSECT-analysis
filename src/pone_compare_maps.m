function T=pone_compare_maps(baseline,alternative,mask,rule)
% Threshold and native-smoothing response summaries.
% rule='threshold' uses the study field-specific zero handling and float32.
% rule='smoothing' accepts already exported maps and their NaNs as validity.
assert(ismember(rule,{'threshold','smoothing'}),'PONE:Rule','Specify threshold or smoothing.');
fields={'DTP_step1','DTP_step2','DTR_step1','DTR_step2','relative_arrival_delay_step1',...
 'steady_state_O2norm_step1','steady_state_O2norm_step2'};
if strcmp(rule,'threshold')
 fields={'refined_lag_step1','relative_arrival_delay_step1','DTP_step1','DTP_step2','DTR_step1','DTR_step2','steady_state_step1','steady_state_step2','steady_state_O2norm_step1','steady_state_O2norm_step2','stimulus_average_step1','stimulus_average_step2','stimulus_average_O2norm_step1','stimulus_average_O2norm_step2'};
end
rows=cell(0,12);
for k=1:numel(fields)
 f=fields{k};a=baseline.(f);b=alternative.(f);
 assert(isequal(size(a),size(b),size(mask)),'PONE:Shape','Aligned maps required.');
 if strcmp(rule,'threshold')
  if ~ismember(f,{'refined_lag_step1','relative_arrival_delay_step1','DTP_step2'}),a(a==0)=NaN;b(b==0)=NaN;end
  a=double(single(a));b=double(single(b));
 end
 va=mask&isfinite(a);vb=mask&isfinite(b);common=va&vb;x=double(a(common));y=double(b(common));d=y-x;
 c=NaN;if numel(x)>1&&std(x)>0&&std(y)>0,q=corrcoef(x,y);c=q(1,2);end
 magnitude=NaN;
 if (startsWith(f,'steady_state')||startsWith(f,'stimulus_average'))&&~isempty(x)
  if strcmp(rule,'threshold')&&median(abs(x))>0,magnitude=100*(median(abs(y))-median(abs(x)))/median(abs(x));end
  if strcmp(rule,'smoothing')&&median(x)~=0,magnitude=100*(abs(median(y))-abs(median(x)))/abs(median(x));end
 end
 rows(end+1,:)={f,nnz(va),nnz(vb),nnz(common),c,median(d),median(abs(d)),sqrt(mean(d.^2)),...
  median(x),median(y),magnitude,nnz(d~=0)};
end
T=cell2table(rows,'VariableNames',{'quantity','baseline_valid','alternative_valid','common_valid',...
 'spatial_r','median_signed_difference','median_absolute_difference','RMSE',...
 'baseline_common_median','alternative_common_median','magnitude_change_percent','changed_voxels'});
end
