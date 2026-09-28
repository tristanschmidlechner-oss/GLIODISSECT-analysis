function R=example_synthetic(outputDir)
% Synthetic normative, patient-deviation and threshold/audit examples.
assert(~isfolder(outputDir)&&~isfile(outputDir),'PONE:OutputExists','Use a new output directory.');mkdir(outputDir);
rng(25851,'twister');G=linspace(-.09,-.03,20)';Y=single((1:40)'/10000+(1:40)'/30*G'+.005*randn(40,20));A=true(size(Y));Y(1,18:20)=NaN;
R=pone_loocv(Y,A,G);pone_write_loocv(R,fullfile(outputDir,'loocv'));
gm=pone_fit_g(Y,A,G);cm=pone_fit_conventional(Y,isfinite(Y));
observed=double(Y(:,2))+.003;mask=true(40,1);
g_deviation=pone_apply_g(observed,G(2),mask,gm);conventional_deviation=pone_apply_conventional(observed,mask,cm);
save(fullfile(outputDir,'synthetic_models.mat'),'gm','cm','g_deviation','conventional_deviation');
[input,mask]=pone_synthetic_input();pairs=[.1 .9;.05 .95;.2 .8];
for k=1:3
 threshold=pone_thresholds(input,mask,pairs(k,:));
 save(fullfile(outputDir,sprintf('threshold_%d.mat',k)),'threshold');
 writetable(threshold.convergence,fullfile(outputDir,sprintf('convergence_%d.csv',k)));
 if k==1,baseline=threshold.maps;else,T=pone_compare_maps(baseline,threshold.maps,mask,'threshold');writetable(T,fullfile(outputDir,sprintf('sensitivity_%d.csv',k)));end
end
end
