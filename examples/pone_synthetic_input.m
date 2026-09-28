function [input,mask]=pone_synthetic_input(seed)
% No participant-derived arrays or fitted parameters.
if nargin<1,seed=25851;end
rng(seed,'twister');N=160;t=1:N;B=60/1.8;W=B;gap=40/1.8;
oxygen=100-50*pulse(t,0,0,B,W,gap);data=zeros(5,2,2,N);
lags=repmat(reshape([0 2 5 8 12],5,1,1),1,2,2);
for k=1:20
 [x,y,z]=ind2sub([5 2 2],k);data(x,y,z,:)=100-3*pulse(t,lags(k),3+mod(k,3)*3,B,W,gap)+.02*randn(1,N);
end
mask=true(5,2,2);input=struct('bBOLD',data,'brainLag',lags,'durationRA',[B W gap W N-B-2*W-gap],...
 'initialGlobalLag',0,'rowScans',t,'shiftedO2',oxygen);
end
function p=pulse(t,tau,rise,B,W,G)
p=ramp(t-B-tau,rise)-ramp(t-B-W-tau,rise)+ramp(t-B-W-G-tau,rise)-ramp(t-B-2*W-G-tau,rise);
end
function y=ramp(x,d)
if d==0,y=double(x>0);else,u=min(max(x/d,0),1);y=(1-cos(pi*u))/2;end
end
