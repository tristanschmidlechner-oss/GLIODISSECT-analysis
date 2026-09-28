function model=pone_fit_conventional(Y,valid,cfg)
% Conventional mean/sample-SD reference with response-map support.
if nargin<3,cfg=pone_config();end
cfg=pone_validate_config(cfg,size(Y,2));
validateattributes(Y,{'single','double'},{'2d','real','nonempty'});
assert(islogical(valid)&&isequal(size(valid),size(Y)),'PONE:Shape','Supply logical map-specific support.');
Y=single(Y);assert(~any(isinf(Y(:))),'PONE:InfiniteInput','Missing values must be NaN.');Y(~valid)=NaN;
n=sum(isfinite(Y),2);mu=mean(Y,2,'omitnan');sd=std(Y,0,2,'omitnan');
eligible=n>=cfg.full_required;mu(~eligible)=NaN;sd(~eligible)=NaN;
p=double(sd(eligible&isfinite(sd)&sd>0));assert(~isempty(p),'PONE:NoVariance','No positive sample SD.');
floorSD=prctile(p,5);sdf=sd;sdf(eligible)=max(sdf(eligible),floorSD);
model=struct('kind','conventional mean/sample SD','mean',mu,'sd_raw',sd,...
 'sd_floored',sdf,'sd_floor',floorSD,'eligible_mask',eligible,'coverage',n,'configuration',cfg);
end
