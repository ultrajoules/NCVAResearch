
%plottop=1;

h1 = figure(30);clf
plot(ee,'o');
%savefig(h1, EigenFigFile)

for plottop=0:0
h2 = figure(31+plottop);clf
fs=9;

nc=3*(1+2*plottop)+2*plottop;
nr=3-plottop*2;

cropxx=0;
cropx=1+cropxx:N-cropxx;

jjj=0;
iispan=[1:3:nr*nc-2*nr*plottop];
iiacumm=0;

for ii=iispan
 jjj=jjj+1;

 %if(ii==iispan(end-2)) jj=29; end;
 %if(ii==iispan(end-1)) jj=13; end;
 %if(ii==iispan(end-0)) jj=15; end;

 if(exist('neigs')&neigs>0)
  eigenvec=(vv(cropx,bb(neigs-(jjj-1))));
 else
  eigenvec=(vv(cropx,bb(2*N-(jjj-1))));
 end
 [temp,mx]=(max(abs(eigenvec')));
 phase0=phase(eigenvec(mx));
 eigenvec=eigenvec/max(max(sqrt(abs(eigenvec).^2)));
 eigenvec=eigenvec*exp(-sqrt(-1)*phase0);
 eigenvecr=real(eigenvec);
 eigenveci=imag(eigenvec);
 eigenveca=abs(eigenvec);


 if(plottop&(mod(jjj,3)~=1))
  iiacumm=iiacumm+1;
 end
 
 eval = ee(bb(end-(jjj-1)));

 subplot(nr,nc,ii+iiacumm)
 plot(y,eigenvecr);
 %title(['$lambda=$',num2str(eval)])
 yl=ylabel({['$\lambda=$',num2str(eval,3)];'real'});
 set(yl,'FontSize',fs)
 set(gca,'FontSize',fs)
 axis tight;ax=axis;axis([ax(1:2) -1 1]);

 subplot(nr,nc,ii+1+iiacumm)
 plot(y,eigenveci);
 ylabel('imag ','FontSize',fs)
 set(gca,'FontSize',fs)
 axis tight;ax=axis;axis([ax(1:2) -1 1]);

 subplot(nr,nc,ii+2+iiacumm)
 plot(y,eigenveca);
 ylabel('abs ','FontSize',fs)
 set(gca,'FontSize',fs)
 axis tight;ax=axis;axis([ax(1:2) 0 1]);

 eval

end
%savefig(h2,EigenVecFigFile)
end %plottop
