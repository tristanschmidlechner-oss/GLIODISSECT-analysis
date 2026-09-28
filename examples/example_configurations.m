function c=example_configurations()
c.study=pone_config();
c.other_cohort=struct('controls',12,'full_required',10,'training_required',9,'q2_required',10);
c.threshold_pairs=[.1 .9;.05 .95;.2 .8];
c.smoothing=struct('spatial_mm',{6,6,4,8},'temporal_volumes',{8,16,12,12},'input_is_unsmoothed',{false,false,true,true});
end
