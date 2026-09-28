function result=pone_apply_conventional(observed,mask,model)
% Floored SD must already be in the patient's grid. No second floor after warp.
assert(strcmp(model.kind,'conventional mean/sample SD'),'PONE:Model','Requires conventional model.');
result=pone_deviation(observed,mask,model.eligible_mask,double(model.mean),double(model.sd_floored),NaN,false);
result.model_definition=model.kind;
end
