function test_pone_loocv()
% Independent backslash-OLS oracle, boundary masks, leakage and output tests.
rng(25851,'twister');G=linspace(-2,2,20)';
Y=single((1:40)'/100 + (1:40)'/30*G' + .15*randn(40,20));
A=true(size(Y));A(1,17:20)=false;A(2,18:20)=false;
Y(3,17:20)=NaN;Y(4,18:20)=NaN;Y(5,:)=single(2);
A(6,:)=false;
R=pone_loocv(Y,A,G);
assert(height(R.folds)==20 && all(R.folds.training_subjects==19));
assert(all(R.folds.training_required==16));
assert(R.voxelwise.heldout_count(1)==0 && R.voxelwise.heldout_count(2)==17);
assert(R.voxelwise.heldout_count(3)==0 && R.voxelwise.heldout_count(4)==17);
assert(all(isnan(R.voxelwise.Q2([1 3 5 6]))));
assert(~isfield(R.in_sample,'Q2_median') && ~isfield(R.cross_validated,'R2_median'));

% Oracle uses per-voxel design matrices, independently of sufficient statistics.
eg=nan(size(Y));em=eg;obs=eg;floorValues=zeros(20,1);
for h=1:20
 train=setdiff(1:20,h);rs=nan(40,1);pg=rs;pm=rs;eligible=false(40,1);
 for v=1:40
  keep=A(v,train)&isfinite(Y(v,train));idx=train(keep);
  if sum(A(v,train))<16 || numel(idx)<16,continue;end
  x=[ones(numel(idx),1),G(idx)];y=double(Y(v,idx))';
  if rank(x)<2,continue;end
  b=x\y;rs(v)=sqrt(sum((y-x*b).^2)/(numel(y)-2));
  % A constant response has exactly zero residual SD analytically; remove
  % backslash roundoff so it is not a spurious positive percentile member.
  if all(y==y(1)),rs(v)=0;end
  eligible(v)=true;pg(v)=[1 G(h)]*b;pm(v)=mean(y);
 end
 floorValues(h)=prctile(rs(eligible&isfinite(rs)&rs>0),5);
 scored=eligible&A(:,h)&isfinite(Y(:,h));
 obs(scored,h)=double(Y(scored,h));
 eg(scored,h)=obs(scored,h)-pg(scored);em(scored,h)=obs(scored,h)-pm(scored);
 near(R.folds.RMSE_G(h),sqrt(mean(eg(scored,h).^2)),1e-10);
 near(R.folds.RMSE_mean(h),sqrt(mean(em(scored,h).^2)),1e-10);
 assert(R.folds.heldout_voxels(h)==nnz(scored));
end
near(R.folds.training_residual_floor,floorValues,1e-9);
for v=1:40
 keep=isfinite(obs(v,:));n=nnz(keep);assert(R.voxelwise.heldout_count(v)==n);
 tss=sum((obs(v,keep)-mean(obs(v,keep))).^2);
 if n>=17 && tss>0
  near(R.voxelwise.Q2(v),1-sum(eg(v,keep).^2)/tss,1e-10);
  near(R.voxelwise.Q2_mean_only(v),1-sum(em(v,keep).^2)/tss,1e-10);
 end
end
pg=sqrt(mean(eg(isfinite(eg)).^2));pm=sqrt(mean(em(isfinite(em)).^2));
near(R.cross_validated.RMSE_G,pg,1e-10);
near(R.cross_validated.RMSE_reduction_pct,100*(pm-pg)/pm,1e-10);
assert(R.cross_validated.pooled_predictions==nnz(isfinite(eg)));
assert(isequal(R.coverage.scored_voxels,R.folds.heldout_voxels));
assert(all(R.coverage.eligible_training_voxels== ...
 R.coverage.scored_voxels+R.coverage.training_eligible_heldout_invalid));

% Arbitrary changes to the held-out response, anatomy and G cannot affect
% that fold's training eligibility or residual-SD percentile.
Y2=Y;Y2(:,1)=single(1000+(1:40)');A2=A;A2(1:10,1)=~A2(1:10,1);
G2=G;G2(1)=100;
R2=pone_loocv(Y2,A2,G2);
assert(R2.folds.training_residual_floor(1)==R.folds.training_residual_floor(1));
assert(R2.folds.eligible_training_voxels(1)==R.folds.eligible_training_voxels(1));
assert(R2.folds.RMSE_G(1)~=R.folds.RMSE_G(1));
% Singular training design must fail closed, never silently fit a slope.
mustError(@()pone_loocv(Y,A,ones(20,1)),'PONE:NoResidualVariance');
mustError(@()pone_loocv(Y(:,1:19),A(:,1:19),G(1:19)),'PONE:ControlCount');
mustError(@()pone_loocv(Y,A(:,1),G),'PONE:Anatomy');
bad=Y;bad(1)=Inf;mustError(@()pone_loocv(bad,A,G),'PONE:InfiniteInput');

% CSV/JSON serialization, separate statistics, no overwrites.
out=tempname;cleanup=onCleanup(@()rmdir(out,'s'));
pone_write_loocv(R,out);
cv=readtable(fullfile(out,'LOOCV_RESULTS.csv'));
ins=readtable(fullfile(out,'NORMATIVE_MODEL_RESULTS.csv'));
assert(ismember('Q2_median',cv.Properties.VariableNames));
assert(~ismember('R2_median',cv.Properties.VariableNames));
assert(ismember('R2_median',ins.Properties.VariableNames));
summary=jsondecode(fileread(fullfile(out,'summary.json')));
near(summary.cross_validated.Q2_median,R.cross_validated.Q2_median,1e-12);
mat=load(fullfile(out,'voxelwise_metrics.mat'));assert(isequaln(mat.voxelwise,R.voxelwise));
assert(isequal(summary.configuration,R.configuration));
mustError(@()pone_write_loocv(R,out),'PONE:OutputExists');
end
function near(a,b,tol)
assert(isequal(size(a),size(b)) && all(abs(a(:)-b(:))<tol),'Independent numeric oracle mismatch.');
end
function mustError(f,id)
try,f();catch e,assert(strcmp(e.identifier,id),'Unexpected error: %s',e.identifier);return;end
error('Expected error was not raised: %s',id);
end
