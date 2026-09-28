function test_pone_models()
rng(12);G=linspace(-2,2,20)';Y=single(.2*randn(30,20)+(1:30)'*G'/10);A=true(size(Y));A(1,16:20)=false;Y(2,17:20)=NaN;
gm=pone_fit_g(Y,A,G);cv=pone_loocv(Y,A,G);
assert(isequaln(gm.in_sample_R2,cv.voxelwise.in_sample_R2));
for v=3:30
 X=[ones(20,1),G];b=X\double(Y(v,:))';assert(max(abs([gm.beta0(v);gm.beta1(v)]-b))<1e-10);
end
cm=pone_fit_conventional(Y,A);V=Y;V(~A)=NaN;
assert(isequaln(cm.mean(3:end),mean(V(3:end,:),2,'omitnan')));
assert(isequaln(cm.sd_raw(3:end),std(V(3:end,:),0,2,'omitnan')));
assert(~cm.eligible_mask(1)&&cm.eligible_mask(2));
obs=double(Y(:,1));m=true(30,1);cg=pone_apply_g(obs,G(1),m,gm);cc=pone_apply_conventional(obs,m,cm);
assert(~isequaln(cc.z,cg.z)&&all(isnan(cc.z(~cc.valid_mask))));
cm.sd_floored(3)=0;cm.sd_floored(4)=-1;cm.sd_floored(5)=single(1e-9);
cc=pone_apply_conventional(obs,m,cm);assert(~cc.valid_mask(3)&&~cc.valid_mask(4)&&cc.valid_mask(5));
gm.residual_sd_raw(3)=-1;cg=pone_apply_g(obs,G(1),m,gm);
assert(abs(cg.z(3)-(obs(3)-gm.beta0(3)-gm.beta1(3)*G(1))/gm.residual_sd_floor)<1e-10);
mustFail(@()pone_apply_conventional(obs,m,gm));mustFail(@()pone_apply_g(obs,NaN,m,gm));
mustFail(@()pone_fit_g(Y,A,ones(20,1)));mustFail(@()pone_fit_conventional(ones(30,20),A));
c=struct('controls',12,'full_required',10,'training_required',9,'q2_required',10);
small=pone_loocv(Y(:,1:12),A(:,1:12),G(1:12),c);assert(height(small.folds)==12&&all(small.folds.training_subjects==11));
% Independent pooled weighting for a non-study cohort with different supports.
counts=small.folds.heldout_voxels;expected=sqrt(sum(counts.*small.folds.RMSE_G.^2)/sum(counts));
assert(abs(expected-small.cross_validated.RMSE_G)<1e-12);
c.training_required=12;mustFail(@()pone_loocv(Y(:,1:12),A(:,1:12),G(1:12),c));
assert(pone_g_scalar([1 NaN -2 0],[true true false true])==.5);
end
function mustFail(f)
try,f();catch,return;end;error('Expected input rejection.');
end
