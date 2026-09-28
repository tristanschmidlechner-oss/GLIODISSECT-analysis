function r=pone_deviation(obs,mask,support,expected,sd,floorValue,isG)
validateattributes(obs,{'single','double'},{'real','nonempty'});
assert(islogical(mask)&&islogical(support),'PONE:Mask','Masks must be logical.');
assert(isequal(size(obs),size(mask),size(support),size(expected),size(sd)),'PONE:Shape','All arrays must be aligned and the same size.');
base=mask&support&isfinite(obs)&isfinite(expected);valid=base&isfinite(sd);
if isG,den=max(sd,floorValue);else,valid=valid&sd>0;den=sd;end
expected(~valid)=NaN;difference=double(obs)-expected;z=difference./den;z(~valid)=NaN;
signed=nan(size(z));signed(valid)=0;signed(valid&z>2)=1;signed(valid&z< -2)=-1;
r=struct('expected',expected,'difference',difference,'z',z,'signed_abnormal',signed,...
 'valid_mask',valid,'nonpositive_input_sd_count',nnz(base&isfinite(sd)&sd<=0));
end
