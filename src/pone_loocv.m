function R = pone_loocv(Y, anatomy, G, cfg)
%PONE_LOOCV Leave-one-out validation of the G-adjusted model.
% Y: voxels x 20, cast to single as in the study map loader.
% anatomy: per-control logical tissue eligibility (GM+WM >= 0.8), same size.
% G: 20 independently calculated native-T1 median response scalars.
% No full-cohort mask/model/floor is accepted or used by a training fold.
% In-sample R2 and cross-validated Q2 are returned in separate structures.
validateattributes(Y,{'single','double'},{'2d','real','nonempty'});
if nargin<4, cfg=pone_config(); end
cfg=pone_validate_config(cfg,size(Y,2));N=cfg.controls;
assert(size(Y,2)==N,'PONE:ControlCount','Control count mismatch.');
assert(islogical(anatomy)&&isequal(size(Y),size(anatomy)), ...
    'PONE:Anatomy','Supply per-control logical anatomy, voxels x 20.');
G=double(G(:));validateattributes(G,{'numeric'},{'numel',N,'finite','real'});
assert(~any(isinf(Y(:))),'PONE:InfiniteInput','Use NaN for missing data, not Inf.');
Y=single(Y);assert(~any(isinf(Y(:))),'PONE:InfiniteInput','Input exceeds single precision range.');nv=size(Y,1);coverage=zeros(N,6);
  valid=anatomy&isfinite(Y);[b0,b1,rs,r2,modelEligible]=pone_fit_training(Y,valid,G,sum(anatomy,2)>=cfg.full_required,cfg.full_required);
  assert(any(modelEligible&isfinite(rs)&rs>0),'PONE:NoResidualVariance','Full model has no positive residual SD.');
  inSampleR2=r2;
  floorModel=prctile(rs(modelEligible&isfinite(rs)&rs>0),5);
  folds=cell(0,15);press=zeros(nv,1);pressmean=press;obsSum=press;obsSq=press;predCount=press;poolN=0;poolS=0;poolM=0;poolBias=0;poolZ=0;poolZn=0;
  for h=1:N
   train=setdiff(1:N,h);[b0t,b1t,rst,~,et,mut]=pone_fit_training(Y(:,train),valid(:,train),G(train),sum(anatomy(:,train),2)>=cfg.training_required,cfg.training_required);
   assert(numel(train)==N-1 && ~ismember(h,train),'PONE:Leakage','Held-out control in training.');
   assert(any(et&isfinite(rst)&rst>0),'PONE:NoResidualVariance','Fold has no positive training residual SD.');
   ft=prctile(rst(et&isfinite(rst)&rst>0),5);ev=et&valid(:,h);ix=find(ev);ob=double(Y(ev,h));pg=b0t(ev)+b1t(ev)*G(h);pm=mut(ev);er=ob-pg;em=ob-pm;zz=er./max(rst(ev),ft);
   assert(all(isfinite(pg))&&all(isfinite(pm))&&all(isfinite(zz)),'PONE:InvalidPrediction','Nonfinite eligible prediction or standardized residual.');
   coverage(h,:)=[nnz(sum(anatomy(:,train),2)>=cfg.training_required), nnz(sum(valid(:,train),2)>=cfg.training_required), nnz(et), nnz(valid(:,h)), nnz(ev), nnz(et&~valid(:,h))];
   press(ix)=press(ix)+er.^2;pressmean(ix)=pressmean(ix)+em.^2;obsSum(ix)=obsSum(ix)+ob;obsSq(ix)=obsSq(ix)+ob.^2;predCount(ix)=predCount(ix)+1;
   poolN=poolN+numel(er);poolS=poolS+sum(er.^2);poolM=poolM+sum(em.^2);poolBias=poolBias+sum(er);poolZ=poolZ+nnz(abs(zz)>2);poolZn=poolZn+nnz(isfinite(zz));
   rg=sqrt(mean(er.^2));rm=sqrt(mean(em.^2));
   folds(end+1,:)={h,nnz(ev),rg,rm,100*(rm-rg)/rm,mean(er),mean(abs(er)),cor(ob,pg),cor(ob,pm),ft,100*mean(abs(zz)>2),N-1,cfg.training_required,nnz(et),rg<rm};
  end
  tss=obsSq-obsSum.^2./max(predCount,1);qe=predCount>=cfg.q2_required&tss>0;q=nan(nv,1);qm=q;q(qe)=1-press(qe)./tss(qe);qm(qe)=1-pressmean(qe)./tss(qe);
  rv=r2(modelEligible&isfinite(r2));qv=q(qe);rmseg=sqrt(poolS/poolN);rmsem=sqrt(poolM/poolN);
  result=struct('model_voxels',nnz(modelEligible),'R2_median',median(rv),'R2_p25',prctile(rv,25),'R2_p75',prctile(rv,75),'residual_floor',floorModel,'Q2_voxels',nnz(qe),'Q2_median',median(qv),'Q2_p25',prctile(qv,25),'Q2_p75',prctile(qv,75),'mean_only_Q2_median',median(qm(qe)),'Q2_positive_pct',100*mean(qv>0),'RMSE_G',rmseg,'RMSE_mean',rmsem,'RMSE_reduction_pct',100*(rmsem-rmseg)/rmsem,'bias',poolBias/poolN,'abs_z_gt2_pct',100*poolZ/poolZn,'pooled_predictions',poolN);
R.configuration=cfg;
R.in_sample=struct('model_voxels',result.model_voxels,'R2_median',result.R2_median, ...
 'R2_p25',result.R2_p25,'R2_p75',result.R2_p75,'residual_floor',result.residual_floor);
R.cross_validated=rmfield(result,{'model_voxels','R2_median','R2_p25','R2_p75','residual_floor'});
R.folds=cell2table(folds,'VariableNames',{'fold','heldout_voxels','RMSE_G','RMSE_mean', ...
 'RMSE_reduction_pct','bias','MAE','spatial_r_G','spatial_r_mean','training_residual_floor', ...
 'abs_z_gt2_pct','training_subjects','training_required','eligible_training_voxels','improved'});
R.cross_validated.controls_improved=sum(R.folds.improved);
R.cross_validated.controls=N;
R.coverage=array2table([(1:N)',coverage],'VariableNames',{'fold','anatomical_training_voxels', ...
 'finite_training_voxels','eligible_training_voxels','valid_heldout_voxels', ...
 'scored_voxels','training_eligible_heldout_invalid'});
R.voxelwise=struct('Q2',q,'Q2_mean_only',qm,'heldout_count',predCount, ...
 'Q2_eligible',qe,'PRESS_G',press,'PRESS_mean',pressmean,'heldout_TSS',tss, ...
 'in_sample_R2',inSampleR2,'in_sample_eligible',modelEligible);
end
function r=cor(a,b)
if numel(a)<3||std(a)==0||std(b)==0;r=NaN;else;c=corrcoef(a,b);r=c(1,2);end
end
