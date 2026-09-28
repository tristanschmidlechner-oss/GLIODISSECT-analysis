function reports=example_synthetic_smoothing(outputDir)
% Complete synthetic NIfTI + oxygen/events example, all supported smoothing configurations.
assert(~isfolder(outputDir)&&~isfile(outputDir),'PONE:OutputExists','Use a new output directory.');mkdir(outputDir);
[input,mask]=pone_synthetic_input(91);paths=cell(160,1);
h=struct('dim',size(mask),'dt',[16 0],'mat',diag([2 2 2 1]),'pinfo',[1;0;0],'descrip','Fully synthetic input');
for k=1:160,h.fname=fullfile(outputDir,sprintf('raw_%04d.nii',k));spm_write_vol(h,input.bBOLD(:,:,:,k));paths{k}=h.fname;end
t=(-5:.5:305)';oxygen=100-50*(double(t>60)-double(t>120)+double(t>160)-double(t>220));
record=zeros(numel(t),23);record(:,1)=1000*t;record(:,6)=oxygen;
endTidal=fullfile(outputDir,'synthetic_EndTidal.txt');writematrix(record,endTidal,'Delimiter','tab');
events=fullfile(outputDir,'synthetic_Events.xlsx');writematrix([0 0;(1:160)'*1800,(1:160)';300000 161],events);
basePaths=pone_spatial_smooth(paths,6,fullfile(outputDir,'baseline_spatial'));
c=example_configurations();configs=[struct('spatial_mm',6,'temporal_volumes',12,'input_is_unsmoothed',false),c.smoothing];
reports=struct([]);
for k=1:numel(configs)
 cfg=configs(k);selected=basePaths;if cfg.input_is_unsmoothed,selected=paths;end
 out=fullfile(outputDir,sprintf('configuration_%d',k));actual=pone_smoothing_workflow(selected,mask,endTidal,events,cfg,out);
 assert(height(actual.convergence)==4&&isfile(fullfile(out,'sensitivity.mat')));
 if cfg.input_is_unsmoothed
  selected=arrayfun(@(n)fullfile(out,'smoothed',sprintf('volume_%04d.nii',n)),1:160,'UniformOutput',false);
 end
 raw=spm_read_vols(spm_vol(char(selected)));raw(raw==0)=NaN;
 filtered=filter_bold_timeseries(raw,mask,1.8,cfg.temporal_volumes);
 timing=prepare_oxygen_timing(endTidal,events,filtered,1.8);lag=estimate_voxel_lag(timing.interpolatedO2,filtered);
 expected=run_hypoxia_bold_core(filtered,lag,timing.durationRA,timing.initialGlobalLag,1:160,@(v)resample_corrected_oxygen(timing,v),mask);
 f=fieldnames(expected);for j=1:numel(f),assert(isequaln(expected.(f{j}),actual.maps.(f{j})),'PONE:SmoothingEquivalence','Driver mismatch %s',f{j});end
 q=struct('spatial_mm',cfg.spatial_mm,'temporal_volumes',cfg.temporal_volumes,'fields_equal',numel(f));
 if isempty(reports),reports=q;else,reports(end+1)=q;end
end
end
