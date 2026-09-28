function result=pone_apply_g(observed,patientG,mask,model)
% Caller supplies aligned coefficients/raw residual SD; floor is applied here.
assert(strcmp(model.kind,'G-adjusted OLS'),'PONE:Model','Requires G-adjusted model.');
validateattributes(patientG,{'numeric'},{'scalar','real','finite'});
validateattributes(model.residual_sd_floor,{'numeric'},{'scalar','finite','positive'});
expected=double(model.beta0)+double(model.beta1)*patientG;
result=pone_deviation(observed,mask,model.eligible_mask,expected,...
 double(model.residual_sd_raw),model.residual_sd_floor,true);
result.model_definition=model.kind;
end
