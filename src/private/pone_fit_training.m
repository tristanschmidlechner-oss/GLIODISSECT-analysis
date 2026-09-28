function [b0,b1,rs,r2,e,mu]=pone_fit_training(Y,V,G,A,minN)
Y=double(Y);Y(~V)=0;n=sum(V,2);sx=double(V)*G;sy=sum(Y,2);sxx=double(V)*(G.^2);sxy=Y*G;syy=sum(Y.^2,2);
cxx=sxx-sx.^2./n;cxy=sxy-sx.*sy./n;sst=syy-sy.^2./n;mu=sy./n;b1=cxy./cxx;b0=mu-b1.*sx./n;sse=max(sst-b1.*cxy,0);
e=A&n>=minN&cxx>0&isfinite(b0)&isfinite(b1);rs=sqrt(sse./(n-2));r2=1-sse./sst;r2(sst<=0)=NaN;
b0(~e)=NaN;b1(~e)=NaN;rs(~e)=NaN;r2(~e)=NaN;mu(~e)=NaN;
end
