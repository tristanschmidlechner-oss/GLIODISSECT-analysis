function model=pone_fit_g(Y,anatomy,G,cfg)
% Voxel x control arrays. Same single-loader/double-OLS arithmetic as LOOCV.
if nargin<4,cfg=pone_config();end
cfg=pone_validate_config(cfg,size(Y,2));
validateattributes(Y,{'single','double'},{'2d','nonempty','real'});
assert(islogical(anatomy)&&isequal(size(anatomy),size(Y)),'PONE:Shape','Logical individual anatomy required.');
G=double(G(:));validateattributes(G,{'double'},{'numel',cfg.controls,'finite','real'});
Y=single(Y);assert(~any(isinf(Y(:))),'PONE:InfiniteInput','Missing values must be NaN.');
V=anatomy&isfinite(Y);
[b0,b1,sd,r2,e]=pone_fit_training(Y,V,G,sum(anatomy,2)>=cfg.full_required,cfg.full_required);
p=sd(e&isfinite(sd)&sd>0);assert(~isempty(p),'PONE:NoVariance','No positive residual SD.');
model=struct('kind','G-adjusted OLS','beta0',b0,'beta1',b1,'residual_sd_raw',sd,...
 'residual_sd_floor',prctile(p,5),'eligible_mask',e,'in_sample_R2',r2,'configuration',cfg);
end
