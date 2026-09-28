function R=pone_smoothing_workflow(inputPaths,mask,endTidalFile,eventsFile,cfg,outputDir)
% User-data workflow. Spatial alternatives require unsmoothed volumes;
% temporal alternatives require the existing 6-mm baseline-smoothed volumes.
assert(isstruct(cfg)&&all(isfield(cfg,{'spatial_mm','temporal_volumes','input_is_unsmoothed'})),'PONE:Config','Invalid smoothing configuration.');
assert(ismember(cfg.spatial_mm,[4 6 8])&&ismember(cfg.temporal_volumes,[8 12 16]),'PONE:Config','Unsupported smoothing.');
assert(cfg.spatial_mm==6||cfg.temporal_volumes==12,'PONE:Config','Joint spatial and temporal alternatives are not supported.');
assert(numel(inputPaths)==160&&islogical(mask),'PONE:Input','160 ordered NIfTI volumes and a logical mask required.');
assert(~isfolder(outputDir)&&~isfile(outputDir),'PONE:OutputExists','Use a new output directory.');
assert(cfg.input_is_unsmoothed||cfg.spatial_mm==6,'PONE:Input','Spatial alternatives must start from unsmoothed data.');
headers=spm_vol(char(inputPaths));
assert(isequal(size(mask),headers(1).dim),'PONE:Shape','Mask must share the volume grid.');
for k=1:160,assert(isequal(headers(k).dim,headers(1).dim)&&max(abs(headers(k).mat(:)-headers(1).mat(:)))<1e-5,'PONE:Geometry','Volume grids differ.');end
mkdir(outputDir);
if cfg.input_is_unsmoothed,inputPaths=pone_spatial_smooth(inputPaths,cfg.spatial_mm,fullfile(outputDir,'smoothed'));end
raw=spm_read_vols(spm_vol(char(inputPaths)));raw(raw==0)=NaN;
input.bBOLD=filter_bold_timeseries(raw,mask,1.8,cfg.temporal_volumes);
timing=prepare_oxygen_timing(endTidalFile,eventsFile,input.bBOLD,1.8);
input.brainLag=estimate_voxel_lag(timing.interpolatedO2,input.bBOLD);
input.durationRA=timing.durationRA;input.initialGlobalLag=timing.initialGlobalLag;
input.rowScans=1:160;input.shiftedO2=@(lag)resample_corrected_oxygen(timing,lag);
R=pone_thresholds(input,mask,[.1 .9]);R.configuration=cfg;
save(fullfile(outputDir,'sensitivity.mat'),'R','-v7.3');
writetable(R.convergence,fullfile(outputDir,'convergence.csv'));
write_maps_nifti_spm(R.maps,inputPaths{1},fullfile(outputDir,'maps'));
end
