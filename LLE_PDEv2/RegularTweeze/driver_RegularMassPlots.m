lw=2;
fs1=18;
fs2=24;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontname', 'Times')
set(0,'DefaultAxesFontSize', fs1)
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultTextFontSize', fs2)
set(0,'DefaultLineLineWidth',lw)
mygreen = [0 0.7 0];

if(exist('QOut')==0) % only load if not loaded
 load('LLERegularMassPlot.mat')
end 
 QOut=QOut;
 QIn = QIn; 
 %smooth data:
 dx=betatau(3)-betatau(2);
 betaN=linspace(betatau(1),betatau(end),length(betatau)*4-3);
 [Xf,B]=meshgrid(cmax,betatau);
 [XfN,BN]=meshgrid(cmax,betaN);
 QOutN=interp2(Xf,B,QOut,XfN,BN);
 QInN=interp2(Xf,B,QIn,XfN,BN);
 %a=601+120;b=501;HopfU16N(a,b)=0; %This is to have a zero in the plot so that colorbar starts at zero

figure(1), clf
axes('Position',[0.08 0.3 0.88 0.5])
colormap(jet(1024))
surfc(XfN,BN,QInN)
shading interp
set(gca,'linewidth',lw,'fontsize',fs1)
xlabel('$\beta$','Interpreter','LaTex','Fontsize',fs2);
ylabel('$\tau_f$','Interpreter','LaTex','Fontsize',fs2); %'Position',[tm-15 0]);
zlabel('$Q_{in}$','Interpreter','LaTex','Fontsize',fs2);
%axis([betatau(1) betatau(end) cmax(1)  cmax(2)])
box on
grid off
%set(gca,'XTick',[0:50:300])
%text(tm+5,2.5,7,'(b)','Color','w','FontSize',fs1)

set(0,'defaulttextinterpreter','none') % This is necessary otehrwise colorbar does not show up when text is in latex
cb=colorbar;
set(cb,'LineWidth',lw)

orient portrait
% for surfaces it is better to save a raster pic otherwise files get huge:
print -djpeg95 -r300 RegularQInC.jpg 
%the following might not work unless you have the appropriate converters:
!convert RegularQIn.jpg -trim RegularQIn.jpg
!jpeg2ps RegularQIn.jpg > RegularQIn.jpg.ps

figure(2), clf
axes('Position',[0.08 0.3 0.88 0.5])
colormap(jet(1024))
surf(XfN,BN,QOutN)
shading interp
set(gca,'linewidth',lw,'fontsize',fs1)
xlabel('$\beta$','Interpreter','LaTex','Fontsize',fs2);
ylabel('$\tau_f$','Interpreter','LaTex','Fontsize',fs2); %'Position',[tm-15 0]);
zlabel('$Q_{out}$','Interpreter','LaTex','Fontsize',fs2);
%axis([betatau(1) betatau(end) cmax(1)  cmax(2)])
box on
grid off
%set(gca,'XTick',[0:50:300])
%text(tm+5,2.5,7,'(b)','Color','w','FontSize',fs1)

set(0,'defaulttextinterpreter','none') % This is necessary otehrwise colorbar does not show up when text is in latex
cb=colorbar;
set(cb,'LineWidth',lw)

orient portrait
% for surfaces it is better to save a raster pic otherwise files get huge:
print -djpeg95 -r300 RegularQOutC.jpg 
%the following might not work unless you have the appropriate converters:
!convert RegularQOut.jpg -trim RegularQOut.jpg
!jpeg2ps RegularQOut.jpg > RegularQOut.jpg.ps

figure(3), clf
axes('Position',[0.08 0.3 0.88 0.5])
colormap(jet(1024))
surf(XfN,BN,QInN)
shading interp
set(gca,'linewidth',lw,'fontsize',fs1)
xlabel('$\beta$','Interpreter','LaTex','Fontsize',fs2);
ylabel('$\tau_f$','Interpreter','LaTex','Fontsize',fs2); %'Position',[tm-15 0]);
zlabel('$Q_{in}$','Interpreter','LaTex','Fontsize',fs2);
axis([betatau(1) betatau(end) cmax(1)  cmax(end)])
box on
grid off
%set(gca,'XTick',[0:50:300])
%text(tm+5,2.5,7,'(b)','Color','w','FontSize',fs1)

set(0,'defaulttextinterpreter','none') % This is necessary otehrwise colorbar does not show up when text is in latex
cb=colorbar;
set(cb,'LineWidth',lw)

orient portrait
% for surfaces it is better to save a raster pic otherwise files get huge:
print -djpeg95 -r300 RegularQIn.jpg 
%the following might not work unless you have the appropriate converters:


figure(4), clf
axes('Position',[0.08 0.3 0.88 0.5])
colormap(jet(1024))
surf(XfN,BN,QOutN)
shading interp
set(gca,'linewidth',lw,'fontsize',fs1)
xlabel('$\beta$','Interpreter','LaTex','Fontsize',fs2);
ylabel('$\tau_f$','Interpreter','LaTex','Fontsize',fs2); %'Position',[tm-15 0]);
zlabel('$Q_{out}$','Interpreter','LaTex','Fontsize',fs2);
axis([betatau(1) betatau(end) cmax(1)  cmax(end)])
box on
grid off
%set(gca,'XTick',[0:50:300])
%text(tm+5,2.5,7,'(b)','Color','w','FontSize',fs1)

set(0,'defaulttextinterpreter','none') % This is necessary otehrwise colorbar does not show up when text is in latex
cb=colorbar;
set(cb,'LineWidth',lw)

orient portrait
% for surfaces it is better to save a raster pic otherwise files get huge:
print -djpeg95 -r300 RegularQOut.jpg 
%the following might not work unless you have the appropriate converters:
