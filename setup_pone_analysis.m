function setup_pone_analysis(gliodissectDir,spmDir)
% Add only the supplement's public folders; private helpers cannot shadow APIs.
assert(strcmp(strtrim(fileread(fullfile(gliodissectDir,'VERSION'))),'1.0.0'),'PONE:Version','GLIODISSECT v1.0.0 required.');
addpath(gliodissectDir);setup_gliodissect;addpath(spmDir);
assert(exist('spm_vol','file')==2,'PONE:SPM','SPM12 required.');
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'src'));addpath(fullfile(root,'examples'));
required={'run_hypoxia_bold_core','estimate_voxel_lag','filter_bold_timeseries','write_maps_nifti_spm'};
for k=1:numel(required)
 p=which(required{k});assert(startsWith(lower(p),lower([gliodissectDir filesep])),'PONE:Path','Dependency function resolved outside the selected release.');
end
end
