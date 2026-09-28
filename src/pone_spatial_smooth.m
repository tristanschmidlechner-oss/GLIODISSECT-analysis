function paths=pone_spatial_smooth(inputPaths,fwhm,outputDir)
% SPM smoothing: unsmoothed, aligned volumes -> float32 NIfTI.
assert(ismember(fwhm,[4 6 8]),'PONE:FWHM','Supported widths are 4, 6, 8 mm.');
assert(iscell(inputPaths)&&~isempty(inputPaths),'PONE:Input','Supply ordered 3-D NIfTI paths.');
assert(~isfolder(outputDir)&&~isfile(outputDir),'PONE:OutputExists','Use a new output directory.');
for k=1:numel(inputPaths),h=spm_vol(inputPaths{k});assert(numel(h)==1,'PONE:Input','One volume per path.');end
mkdir(outputDir);paths=cell(size(inputPaths));
for k=1:numel(inputPaths)
 paths{k}=fullfile(outputDir,sprintf('volume_%04d.nii',k));
 spm_smooth(inputPaths{k},paths{k},[fwhm fwhm fwhm],16);
end
end
