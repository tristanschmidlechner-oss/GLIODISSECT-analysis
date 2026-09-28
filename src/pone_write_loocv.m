function pone_write_loocv(R, outputDir)
%PONE_WRITE_LOOCV Machine-readable results; caller controls private output.
% Fold ordinals contain no source identifiers. Never put real-data voxel outputs
% or subject-level tables in a public release without a separate disclosure review.
assert(~isfolder(outputDir)&&~isfile(outputDir),'PONE:OutputExists','Use a new output directory.');
mkdir(outputDir);
writetable(struct2table(R.cross_validated),fullfile(outputDir,'LOOCV_RESULTS.csv'));
writetable(struct2table(R.in_sample),fullfile(outputDir,'NORMATIVE_MODEL_RESULTS.csv'));
writetable(R.folds,fullfile(outputDir,'LOOCV_subjects.csv'));
writetable(R.coverage,fullfile(outputDir,'LOOCV_coverage.csv'));
summary=struct('cross_validated',R.cross_validated,'in_sample',R.in_sample,'configuration',R.configuration);
fid=fopen(fullfile(outputDir,'summary.json'),'w');assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));fwrite(fid,jsonencode(summary),'char');
voxelwise=R.voxelwise;save(fullfile(outputDir,'voxelwise_metrics.mat'),'voxelwise','-v7.3');
end
