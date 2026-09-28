function R=example_private_arrays(inputMat,newOutputDir,cfg)
% MAT variables: Y (voxels x controls), anatomy (logical, same size), G (controls x 1).
% Supply a configuration explicitly for cohorts other than the 20-control study.
if nargin<3,cfg=pone_config();end
q=load(inputMat,'Y','anatomy','G');R=pone_loocv(q.Y,q.anatomy,q.G,cfg);pone_write_loocv(R,newOutputDir);
end
