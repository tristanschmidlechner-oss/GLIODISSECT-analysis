function test_pone_spm()
root=tempname;mkdir(root);c=onCleanup(@()rmdir(root,'s'));
h=struct('fname',fullfile(root,'synthetic.nii'),'dim',[5 4 3],'dt',[16 0],...
 'mat',diag([2 2 2 1]),'pinfo',[1;0;0],'descrip','Synthetic test');
rng(44);x=100+randn(h.dim);spm_write_vol(h,x);
for width=[4 6 8]
 out=pone_spatial_smooth({h.fname},width,fullfile(root,sprintf('s%d',width)));
 direct=fullfile(root,sprintf('direct%d.nii',width));spm_smooth(h.fname,direct,[width width width],16);
 assert(isequaln(spm_read_vols(spm_vol(out{1})),spm_read_vols(spm_vol(direct))));
end
[in,mask]=pone_synthetic_input();
for span=[8 12 16]
 filtered=filter_bold_timeseries(in.bBOLD,mask,1.8,span);assert(isequal(size(filtered),size(in.bBOLD))&&all(isfinite(filtered(:))));
end
R=pone_thresholds(in,mask);h.dim=size(mask);h.fname=fullfile(root,'reference.nii');spm_write_vol(h,zeros(size(mask)));
out=write_maps_nifti_spm(R.maps,h.fname,fullfile(root,'maps'));assert(numel(out)==13);
v=spm_read_vols(spm_vol(out{12}));assert(isequaln(v,double(single(R.maps.steady_state_O2norm_step1))));
end
