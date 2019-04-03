set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontname', 'Times')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultTextFontSize', 12)
set(0,'DefaultLineLineWidth',1)


global allt allx allu alpha mu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NLSE variable definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global delta T0 eta P ep V 
load('/Users/jules/Documents/Jules/BEC Research /XuCoen/Code/Stability05112015/SSB_Stability_UnStableBranch_Best.mat')
load('SSB_Stability_NEW_02072015.mat')
uiniS = uiniSave;
VV = Vsave;
eeSave=[];
vvSave=[];
uiniSave =[];
%tgrowthSave =[];
%growthSave=[];
Vsave=[];
%uPSave=[];
%paramSave =[];

ep1 = 1; %0.31624; %0; %1; %0.31624*0;
ep = 1; %0.1; %0.1; %0.1; %1; %0.1;
%aa = 0.5;  %add to shift spectrum - such that equation is -i(1-aa)
delta = 0.92;
T0 = 2.3; 
eta = -1;
onehalf = 1.0;  %FOr XuCoen should be 1
g = -1; %Bright/focusing/attractive = -1,   Dark/defocusing/repulsive = +1
alpha = g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params for Newton: 
it=15; % Newton iterations for each fixed point
normconv=1e-95; % Tolerance for convergence of Newton
mtalpha = -3*alpha;

% Spatial computational domain defined
L=15;  % domain
dx = 0.1;
x = [-L:dx:L];
N = length(x);
y=x';
%dx=y(2)-y(1);
onehalf = 1; %0.5;         % NLS: (onehalf)*u_xx or (1)*u_xx
onehalfoodx2=onehalf/dx^2;
xi0 = x';
oodx2=1/dx^2;
oos = 1/6;
dxx = 1/dx^2;
ohdxx = onehalf*dxx;
onehalfoodx2=onehalf/(dx*dx);

%uini = 1e-6 + uiniSave(219,:)';
P = 6.4;
    uini = uiniS(128,:)';
    128
uini = uiniS(94,:)';
        
uini = uiniS(300,:)';

for i = 93 %340:-1:300; %300:-1:210; %128:length(PSave) %128:-1:91;
    vvv =[];
    eee=[];
    i
    P = PSave(i);
    uini = 0.5*uiniS(i,:)';
    %uini = uiniS(128,:)'; 
    %uini = 0.35*uini.*uiniS(i,:)';
    %uini = uiniS(i+1,:)';
    V = VSave(i,:)';
    %S = @(x, P) sqrt(P).*exp(-(x/T0).^2);
    %V = S(xi0, P);
    allAumf = [];
    neigs=100; % # of evals to compute. neigs=0 => FULL spectrum
    allee=[];

    u0 = [reshape(real(uini),N,1); reshape(imag(uini),N,1); 0];

    alpha = -1;
    maxit=1000;                          %nsoli params
    maxitl=1000;
    etamax=0.9;
    lmeth=3;
    restart_limit=1000;
    sol_parms=[maxit,maxitl,etamax,lmeth,restart_limit];
    [uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,eta,dx,xi0),1e-15*[1,1],sol_parms);
    RHS = GPE_nsoli_Stationary_1D(uNSOLI,N,onehalfoodx2,V,eta,dx,xi0);
    figure(500);
    plot(1:2*N+1, RHS)
    
    ubefore = reshape(uNSOLI(1:N) + 1i*uNSOLI([1:N]+N),N,1);
    figure(15);
    plot(xi0, real(ubefore), xi0, imag(ubefore), xi0, ubefore.*conj(ubefore))
    u = ubefore;
    uini = u;
    
    stability1D
    [v0, e0] = eigs(J,1,0);
    [v2, e2] = eigs(J,1,-2);
    
    eeSave(i,:)=ee;
    vvSave(i,:,:)=vv;
    uiniSave(i,:)=uini;
    Vsave(i,:)=V;
    pitchforkStability_0(i,:)=e0;
    pitchforkStability_2(i,:)=e2;
    pitchforkStability_v0(i,:)=v0;
    pitchforkStability_v2(i,:)=v2;
    %growthSave(i,:)=growth;
    %uPSave(i,:,:)=uP;
    %paramSave(i,:) =param;

end

jj=0;
eeUnstable =[];
PUnstable_New =[];
for i = 1:length(eeSave)
    for k = 1:100
if (real(eeSave(i,k)) > 0)
    jj = jj+1;
    eeUnstable(jj) = real(eeSave(i,k));
    PUnstable_New(jj) = PSave(i);
end
    end
end


    eeSaveStable(i,:)=[];
    vvSaveStable(i,:,:)=[];
    uiniSaveStable(i,:)=[];
    VsaveStable(i,:)=[];
    pitchforkStability_0Stable(i,:)=[];
    pitchforkStability_2Stable(i,:)=[];
    pitchforkStability_v0Stable(i,:)=[];
    pitchforkStability_v2Stable(i,:)=[];
PSaveStable = [4.50:0.05:10.65];
load('/Users/jules/Documents/Jules/BEC Research /XuCoen/Code/XuCoenGaussian6/NSOLISolution_SSP6.4.mat')
uini=ubefore;
%39
for i = 39:-1:1; %64:length(PSaveStable) %340:-1:300; %300:-1:210; %128:length(PSave) %128:-1:91;
    vvv =[];
    eee=[];
    i
    P = PSaveStable(i);
    %uini = uiniS(128,:)'; 
    %uini = 0.35*uini.*uiniS(i,:)';
    %uini = uiniS(i+1,:)';
    %V = VV(i,:)';
    S = @(x, P) sqrt(P).*exp(-(x/T0).^2);
    V = S(xi0, P);
    allAumf = [];
    neigs=100; % # of evals to compute. neigs=0 => FULL spectrum
    allee=[];

    u0 = [reshape(real(uini),N,1); reshape(imag(uini),N,1); 0];

    alpha = -1;
    maxit=1000;                          %nsoli params
    maxitl=1000;
    etamax=0.9;
    lmeth=3;
    restart_limit=1000;
    sol_parms=[maxit,maxitl,etamax,lmeth,restart_limit];
    [uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,eta,dx,xi0),1e-15*[1,1],sol_parms);
    RHS = GPE_nsoli_Stationary_1D(uNSOLI,N,onehalfoodx2,V,eta,dx,xi0);
    figure(500);
    plot(1:2*N+1, RHS)
    
    ubefore = reshape(uNSOLI(1:N) + 1i*uNSOLI([1:N]+N),N,1);
    figure(15);
    plot(xi0, real(ubefore), xi0, imag(ubefore), xi0, ubefore.*conj(ubefore))
    u = ubefore;
    uini = u;
    
    stability1D
    [v0, e0] = eigs(J,1,0);
    [v2, e2] = eigs(J,1,-2);
    
    eeSaveStable(i,:)=ee;
    vvSaveStable(i,:,:)=vv;
    uiniSaveStable(i,:)=uini;
    VsaveStable(i,:)=V;
    pitchforkStability_0Stable(i,:)=e0;
    pitchforkStability_2Stable(i,:)=e2;
    pitchforkStability_v0Stable(i,:)=v0;
    pitchforkStability_v2Stable(i,:)=v2;
    %growthSave(i,:)=growth;
    %uPSave(i,:,:)=uP;
    %paramSave(i,:) =param;

end    

save('SSB_Stability_NEW_02072015.mat','eeSave','vvSave','uiniSave','Vsave','eeUnstable','PUnstable_New','PSave','xi0',...
    'pitchforkStability_0', 'pitchforkStability_2', 'pitchforkStability_v0', 'pitchforkStability_v2',...
    'eeSaveStable','vvSaveStable','uiniSaveStable','VsaveStable','pitchforkStability_0Stable','pitchforkStability_2Stable',...
    'pitchforkStability_v0Stable','pitchforkStability_v2Stable')

   
    % pert=1e-6;
    % rand('twister',5489);
    % %u = u+pert*(rand(N,1)-0.5);
    % u = u+pert*vv(:,1);
    %
    % dt=0.001;tf=10;snaps=100;disp=25;
    % t0=0;Vini=V'; damping_step = 0;
    % Nx=N;Lx=L;
    % [uP,uPt] = cNLS1D(snaps,100,tf/dt+1,N,dt,xi0,u,ohdxx,eta,V,oos,allx,allt,damping_step,dx);
    %
    % alltgrowth=[];
    % growth=[];
    %
    % t=0;
    % j=1;
    % for kk=1:length(uPt);
    %     u = uP(kk,:)';
    %     t = uPt(kk);
    %     Mu=max(uini.*conj(uini));
    %     mMu=min(min(uini.*conj(uini)));
    %     h3 = figure(1); clf;
    %     subplot(3,1,1)
    %     set(gca,'FontSize',[16]);
    %     mytext=['T=(' num2str(t0) '+' num2str(kk*dt) ')/(' num2str(t0) '+' num2str(tf) ')'];
    %     plot(x,u.*conj(u),x,V.*conj(V),x,uini.*conj(uini),'--k');
    %     axis([-Lx Lx -0.1 2.0]);
    %     axis([-Lx Lx mMu-0.1 Mu*1.1]);
    %     title(['t=',num2str(t)])
    %     subplot(3,2,3)
    %     plot(x,u.*conj(u)-uini.*conj(uini));
    %     ylabel('pert')
    %     axis tight
    %     subplot(3,2,4)
    %     plot(x,real(u),x,imag(u),x,abs(u));
    %     ylabel('Re,Im,Abs')
    %     axis([-Lx Lx -sqrt(Mu)*1.1 sqrt(Mu)*1.1]);
    %     alltgrowth=[alltgrowth,t];
    %     growth=[growth,sum(abs(u.*conj(u)-uini.*conj(uini)))*dx];
    %     subplot(3,2,5)
    %     plot(alltgrowth,abs(growth-growth(1)),'o-');
    %     xlabel('t');ylabel('|pert|');
    %     axis tight;
    %     subplot(3,2,6)
    %     semilogy(alltgrowth(2:end),abs(growth(2:end)-1*growth(1)),'o-');
    %     xlabel('t');ylabel('|pert|');
    %     axis tight; ax=axis;
    %     if(length(growth)>disp/6)
    %         lg=length(growth);
    %         iii=floor(2*lg/4):lg;
    %         xdata=alltgrowth(iii)-alltgrowth(iii(1));
    %         ydata=log(growth(iii));
    %         param=real(polyfit(xdata,ydata,1));
    %         hold on
    %         semilogy(alltgrowth(iii),exp(param(2))*exp((alltgrowth(iii)-alltgrowth(iii(1)))*param(1)),'r-');
    %         hold off
    %         ax(3)=min(growth(2)-growth(1));
    %         if(ax(4)-ax(3))
    %             axis(ax);
    %         else
    %             axis tight;
    %         end
    %         title(['$\lambda=$',num2str(param(1),8)])
    %     end
    %     drawnow
    % end
    % savefig(h3,PertFigFile);
    
    
