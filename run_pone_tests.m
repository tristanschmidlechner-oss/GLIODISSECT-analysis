function results=run_pone_tests()
% Run after setup_pone_analysis. Only synthetic arrays/files are created.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'tests'));c=onCleanup(@()rmpath(fullfile(root,'tests')));
test_pone_loocv;fprintf('PASS independent LOOCV oracle, isolation, validity and serialization\n');
test_pone_models;fprintf('PASS distinct models, coverage, non-study cohort and invalid inputs\n');
results.adapter=pone_verify_adapter();fprintf('PASS release/study-core adapter equivalence\n');
test_pone_spm;fprintf('PASS SPM smoothing, temporal spans and NIfTI round-trip\n');
out=tempname;c2=onCleanup(@()cleanup(out));results.smoothing=example_synthetic_smoothing(out);
fprintf('PASS complete synthetic smoothing workflow, five configurations\n');
results.test_groups=5;results.passed=true;results.matlab=version;results.spm=spm('Ver');
end

function cleanup(out)
if isfolder(out),rmdir(out,'s');end
end
