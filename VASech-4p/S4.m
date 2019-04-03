clear all; close all; 
global delta h sigma_X Ein tau0
%In general, you can generate N random numbers in the 
%interval [a,b] with the formula r = a + (b-a).*rand(N,1).%
rng(0,'twister');
N = 1000;
%choose a [-5 5]
a_random = 0 + (10).*rand(N,1);
%choose b [-5 5]
b_random = -5 + (5+5).*rand(N,1);
%choose c [-3 3]
c_random =  -5 + (5+5).*rand(N,1);
%choose d [-3 3]
d_random =  0 + (3).*rand(N,1);
%choose sigma [0.1 3]
sigma_random = 0.1 + (3-0.1).*rand(N,1);
%choose sigma [0.1 3]
xi_random = -10 + (20).*rand(N,1);

Ein = sqrt(.960)*sqrt(1.2*0.1*0.1/0.146^3);
h =  2; %2;
sigma_X = 2;
tau0 = 2;
delta =  0.41/0.146; 

L = 20;
dx = 0.1;
x = (-L:dx:L);
xi0 = x';
fitting_func =  inline(...
        '[(p(1)).*sech(p(1)*(x-p(4))).*exp(1i*(p(3).*(x-p(4)) + p(2)))]'...
        ,'p','x');
param0 = [sqrt(3), 1, 1 , 0];
paramSeed=param0; %lsqcurvefit(fitting_func,param0,xi0, real(uPDE),[],[]);
uSeed = (paramSeed(1)).*sech(paramSeed(1)*(xi0-paramSeed(4))).*exp(1i*(paramSeed(3).*(xi0-paramSeed(4)) + paramSeed(2)));
a_random(1) = real(paramSeed(1));
b_random(1) = real(paramSeed(2));
c_random(1) = real(paramSeed(3));
xi_random(1) = real(paramSeed(4));
plot( xi0, uSeed.*conj(uSeed))

maxit=1000;                          %nsoli params
maxitl=1000;
etamax=0.9;
lmeth=3;
restart_limit=1000;
sol_parms=[maxit,maxitl,etamax,lmeth,restart_limit];

% options.Display= 'final-detailed'; %'iter';
% options.LargeScale = 'off';
% options.TrustRegionReflective = 'On';
% options.Hessian = 'On';
% options.DerivativeCheck = 'On';
% options.LineSearchType = 'quadcubic';
options.MaxFunEvals = 10000;
options.MaxIter = 10000;
options.TolFun = 1e-8;
options.TolX = 1e-8;

for i = 2; %:10;
    x0 = [a_random(i); b_random(i); c_random(i); xi_random(i)];
    [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_S4(solp_),1e-15*[1,1],sol_parms);
    %[uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,eta,dx,xi0),1e-15*[1,1],sol_parms);
    %[solp, resnorm,Ein residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[],options);
    %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
    uVA(:,i) = (solp(1)).*sech(solp(1)*(xi0-solp(4))).*exp(1i*(solp(3).*(xi0-solp(4)) + solp(2)));
    %thetauVA = atan2(imag(uVA(:,i)), real(uVA(:,i)));
    %uVA(:,i) = uVA(:,i).*exp(-1i*thetauVA);
    params(:,i) = [solp(1); solp(2); solp(3); solp(4)];
    %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
    bVA(i) = uVA(151,i).*conj(uVA(151,i));
    fRHS(:,i) = LLE_steady_state_S4(solp);
    error(i) = ierr; %exitflag;
end

uVA2 = [];
pVA =[];
test =[];
tmax = 25;
t = 0:0.01:tmax;
p0 = params(:,2);
[VAt,xVA]=ode45('VAF',t,[p0]);

xi0 = xi0;
for i = 1:length(t);
    pVA = xVA(i,:);
    uVA2(:,i) = (pVA(1)).*sech(pVA(1)*(xi0-pVA(4))).*exp(1i*(pVA(3).*(xi0-pVA(4)) + pVA(2)));
    %plot(xi0, uVA2(:,i).*conj(uVA2(:,i)))
    test(:,i) = uVA2(:,i).*conj(uVA2(:,i));
    %drawnow;
end
figure(215)
colormap jet
mesh(t,xi0,test)
shading interp
axis([0 tmax -5 5])

% N  = 788;
% %load('XC6_nSoli_P64_ParamTest.mat')
% figure(50)
% hold on;
% scatter(P, bPDE);
% for i = 1:N
%     scatter(P, bVA(i));
% end
% hold off;
% 
% figure(10)
% plot(xi0, uPDE.*conj(uPDE), xi0, uVA(:,593).*conj(uVA(:,593)),xi0, uVA(:,467).*conj(uVA(:,467)),xi0, uVA(:,281).*conj(uVA(:,281)));
% xlim([-5, 5]);
% legend('PDE', 'VA Best', 'VA')
% 
% figure(150)
% scatter((1:N), leastSquares);
% ylim([0,200]);
% figure(151)
% scatter((1:N), bVA)
% 
% %Q = -iu + iS
% %A = 593;
% % A=1;
% % B= 32;
% % C = 70;
% % D = 32;
% %B = 467; 
% %C = 281;
% %D = 944;
% 
% 
% %load('XC6_delta_P64_ParamTest.mat');
% %Delta  Q = -iu + Delta u + iS
% A = 535;
% B = 440;
% C = 120;
% D = 284;
% D = 552;
% 
% figure(100)
% plot(xi0, uPDE.*conj(uPDE), xi0, uVA(:,A).*conj(uVA(:,A)),xi0, uVA(:,B).*conj(uVA(:,B)),xi0, uVA(:,C).*conj(uVA(:,C)));
% xlim([-5, 5]);
% legend('PDE', 'VA Best', 'VA')
% 
% Pump = (0:0.1:14);
% x0 = params(:,A);
% NN = 65;
% %NN = 20;
% %for j = 20:length(Pump)
% %for j = NN:length(Pump)
% for j = NN:length(Pump)
% 
%     P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B1(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B1(:,j) = [solp(1); solp(2); solp(3); solp(4); solp(5); solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B1(j) = uVA_B1(151,j).*conj(uVA_B1(151,j));
%     fRHS_B1(:,j) = LLE_steady_state_XC6(solp);
%     error_B1(j) = ierr; %exitflag;
%     x0 = params_B1(:,j);
% end
% x0 = params(:,A);
% %for j = 19:-1:1;
% for j = NN-1:-1:1;
%      P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B1(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B1(:,j) = [solp(1); solp(2); solp(3); solp(4);solp(5);solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B1(j) = uVA_B1(151,j).*conj(uVA_B1(151,j));
%     fRHS_B1(:,j) = LLE_steady_state_XC6(solp);
%     error_B1(j) = ierr; %exitflag;
%     x0 = params_B1(:,j);
% end
% 
% 
% x0 = params(:,B);
% 
% %NN = 20;
% %for j = 20:length(Pump)
% for j = NN:length(Pump)
%     P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B2(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B2(:,j) = [solp(1); solp(2); solp(3); solp(4); solp(5); solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B2(j) = uVA_B2(151,j).*conj(uVA_B2(151,j));
%     fRHS_B2(:,j) = LLE_steady_state_XC6(solp);
%     error_B2(j) = ierr; %exitflag;
%     x0 = params_B2(:,j);
% end
% x0 = params(:,B);
% %for j = 19:-1:1;
% for j = NN-1:-1:1;
%      P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B2(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B2(:,j) = [solp(1); solp(2); solp(3); solp(4);solp(5);solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B2(j) = uVA_B2(151,j).*conj(uVA_B2(151,j));
%     fRHS_B2(:,j) = LLE_steady_state_XC6(solp);
%     error_B2(j) = ierr; %exitflag;
%     x0 = params_B2(:,j);
% end
% 
% x0 = params(:,C);
% 
% %NN = 20;
% %for j = 20:length(Pump)
% for j = NN:length(Pump)
%     P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B3(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B3(:,j) = [solp(1); solp(2); solp(3); solp(4); solp(5); solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B3(j) = uVA_B3(151,j).*conj(uVA_B3(151,j));
%     fRHS_B3(:,j) = LLE_steady_state_XC6(solp);
%     error_B3(j) = ierr; %exitflag;
%     x0 = params_B3(:,j);
% end
% x0 = params(:,C);
% %for j = 19:-1:1;
% for j = NN-1:-1:54; %1;
%      P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
% %     [solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B3(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B3(:,j) = [solp(1); solp(2); solp(3); solp(4);solp(5);solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B3(j) = uVA_B3(151,j).*conj(uVA_B3(151,j));
%     fRHS_B3(:,j) = LLE_steady_state_XC6(solp);
%     error_B3(j) = ierr; %exitflag;
%     x0 = params_B3(:,j);
% end
% 
% x0 = params(:,D);
% 
% %NN = 20;
% %for j = 20:length(Pump)
% for j = NN:length(Pump)
%     P = Pump(j);
%     [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B4(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B4(:,j) = [solp(1); solp(2); solp(3); solp(4); solp(5); solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B4(j) = uVA_B4(151,j).*conj(uVA_B4(151,j));
%     fRHS_B4(:,j) = LLE_steady_state_XC6(solp);
%     error_B4(j) = ierr; %exitflag;
%     x0 = 0.1 + params_B4(:,j);
% end
% x0 = params(:,D);
% %for j = 19:-1:1;
% for j = NN-1:-1:1;
%      P = Pump(j);
%      [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_XC6(solp_),1e-15*[1,1],sol_parms);
%     %[solp, resnorm, residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[]);
%     %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
%     uVA_B4(:,j) =  (solp(1)).*exp(-0.5*(((xi0 - solp(6)).^2)/solp(5).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
%     params_B4(:,j) = [solp(1); solp(2); solp(3); solp(4);solp(5);solp(6)];
%     %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
%     bVA_B4(j) = uVA_B4(151,j).*conj(uVA_B4(151,j));
%     fRHS_B4(:,j) = LLE_steady_state_XC6(solp);
%     error_B4(j) = ierr; %exitflag;
%     x0 = 0.1 + params_B4(:,j);
% end
% 
% 
% load('SSBBifurcationBranch.mat')
% sPumpPDE = [PSave(22:131); PSave(163:193)];
% sBranchPDE = [bifUSave(22:131);bifUSave(163:193)]; 
% load('SSBBifurcationBranch4.mat')
% aSPumpPDE = [PSave(45:52);PSave(55:76); PSave(5:8); PSave(84); PSave(10:19); PSave(89); PSave(21:24); PSave(85); PSave(26:31); PSave(94); PSave(33:35)];
% asBranchPDE = [bifUSave(45:52);bifUSave(55:76);bifUSave(5:8); bifUSave(84); bifUSave(10:19); bifUSave(89); bifUSave(21:24); bifUSave(85); bifUSave(26:31); bifUSave(94); bifUSave(33:35)];
% load('../XuCoenGaussian6/XC6G_P64_Branch.mat')
% bVA_B1(1)=0; 
% bVA_B3(1)=0;
% clf; 
% figure(300)
% hold on;
% plot(Pump, bVA_B3,'g',Pump, bVA_B1,'b','Linewidth',2)
% %plot(Pump, bVA_B3,Pump, bVA_B4,'Linewidth',2)
% %Pump, bVA_B2 Pump, bVA_B4, Pump, bVA_B1 ,Pump, bVA_B4,
% %plot(Pump, bVA_B3, Pump, bVA_B4,'Linewidth',2)
% %legend('PDE','PDE Antisymmetric');
% plot(sPumpPDE, sBranchPDE, 'k', aSPumpPDE, asBranchPDE, 'r', 'Linewidth',2)
% %title('\bf{SSB Bifurcation Diagram}','Interpreter','LaTex','Fontsize',36)
% xlabel('$X$','Interpreter','LaTex','Fontsize',24);
% ylabel('$|u(\tau =0)|^2$','Interpreter','LaTex','Fontsize',24);
% %legend('NCVA - Asymmetric','NCVA - Symmetric ','PDE - Symmetric','PDE - Asymmetric');
% 
% load('../PDENSOLiSolutions/NSOLISolution_P4.mat')
% allx = xi0;
% SS4 = sqrt(4)*exp(-(xi0/T0).^2);
% uVA4 = uVA_B1(:,41);
% u4 = u;
% theta = atan2(imag(u4), real(u4));
% u4 = u4.*exp(-1i*theta);
% load('../PDENSOLiSolutions/NSOLISolution_P6.4.mat')
% S64 = sqrt(6.4)*exp(-(xi0/T0).^2);
% uVAS65 = uVA_B1(:,65);
% uS65 = u;
% theta = atan2(imag(uS65), real(uS65));
% uS65 = uS65.*exp(-1i*theta);
% load('../PDENSOLiSolutions/NSOLISolution_Bi2P6.4.mat')
% uVAaS65 =uVA_B3(:,65);
% uaS65 = u;
% theta = atan2(imag(uaS65), real(uaS65));
% uaS65 = uaS65.*exp(-1i*theta);
% load('../PDENSOLiSolutions/NSOLISolution_P10.5.mat')
% S105 = sqrt(10.5)*exp(-(xi0/T0).^2);
% uVA105 = uVA_B1(:,106);
% uVA105as = uVA_B3(:,106);
% u105 = u;
% theta = atan2(imag(u105), real(u105));
% u105 = u105.*exp(-1i*theta);
% %Create the new axes
% % ax1 = axes('Position',get(ax1,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','left',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % ax2 = axes('Position',get(ax2,'Position'),...
% %            'XAxisLocation','top',...
% %            'YAxisLocation','middle',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % ax3 = axes('Position',get(ax3,'Position'),...
% %            'XAxisLocation','bottom',...
% %            'YAxisLocation','middle',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% % ax4 = axes('Position',get(ax4,'Position'),...
% %            'XAxisLocation','bottom',...
% %            'YAxisLocation','right',...
% %            'Color','none',...
% %            'XColor','k','YColor','k');
% %ax1 = axes('position',[.25 .2 0.15 0.15],...
% ax1 = axes('position',[.23 .2 0.15 0.15],...
%             'Color','none',...
%             'XColor','k','YColor','k');
% %plot(allx, u4.*conj(u4), 'k', allx, uVA4.*conj(uVA4),'b',allx, S4.*conj(S4), 'y.','Linewidth',2)
% plot(allx, u4.*conj(u4), 'k', allx, uVA4.*conj(uVA4),'b','Linewidth',2)
% xlabel('Fast-time $\tau$','Interpreter','LaTex','Fontsize',12);
% ylabel('$|u(\tau)|^2$','Interpreter','LaTex','Fontsize',12);
% xlim([-5 5]);
% title('$X = 4$','Interpreter','LaTex','Fontsize',24) 
% ax2 = axes('position',[.35 .74 0.15 0.15],...
%             'Color','none',...
%             'XColor','k','YColor','k');
% %plot(allx, uS65.*conj(uS65), 'k', allx, uVAS65.*conj(uVAS65),'b',allx, S64.*conj(S64), 'y.', 'Linewidth',2)
% plot(allx, uS65.*conj(uS65), 'k', allx, uVAS65.*conj(uVAS65),'b', 'Linewidth',2)
% xlabel('Fast-time $\tau$','Interpreter','LaTex','Fontsize',12);
% ylabel('$|u(\tau)|^2$','Interpreter','LaTex','Fontsize',12);
% xlim([-5 5]);
% title('$X = 6.4$','Interpreter','LaTex','Fontsize',24) 
% ax3 = axes('position',[.57 .32 0.15 0.15],...
%             'Color','k',...
%             'XColor','k','YColor','k');
% %plot(allx, uaS65.*conj(uaS65), 'r', allx, uVAaS65.*conj(uVAaS65),'g',allx, S64.*conj(S64), 'y.', 'Linewidth',2)
% plot(allx, uaS65.*conj(uaS65), 'r', allx, uVAaS65.*conj(uVAaS65),'g', 'Linewidth',2)
% xlabel('Fast-time $\tau$','Interpreter','LaTex','Fontsize',12);
% ylabel('$|u(\tau)|^2$','Interpreter','LaTex','Fontsize',12);
% xlim([-5 5]);
% title('$X = 6.4$','Interpreter','LaTex','Fontsize',24) 
% ax4 = axes('position',[.75 .5 0.15 0.15],...
%             'Color','k',...
%             'XColor','k','YColor','k');
% %plot(allx, u105.*conj(u105), 'k', allx, uVA105.*conj(uVA105),'b',allx, uVA105as.*conj(uVA105as),'g',allx, S105.*conj(S105),'y.','Linewidth',2)
% %plot(allx, u105.*conj(u105), 'k', allx, uVA105.*conj(uVA105),'b',allx, S105.*conj(S105),'y.','Linewidth',2)
% plot(allx, u105.*conj(u105), 'k', allx, uVA105.*conj(uVA105),'b','Linewidth',2)
% xlabel('Fast-time $\tau$','Interpreter','LaTex','Fontsize',12);
% ylabel('$|u(\tau)|^2$','Interpreter','LaTex','Fontsize',12);
% xlim([-5 5]);
% title('$X = 10.5$','Interpreter','LaTex','Fontsize',24) 
% hold off;
% 
% figure(3)
% hold on;
% for i = 1:6;
% plot(Pump, fRHS_B1(i,:));
% end
% hold off;
% 
% save('XC6G_P64_ParamTest.mat', 'uVA', 'params', 'leastSquares', 'bVA', 'fRHS','error'); 
% save('XC6G_P64_Branch.mat','Pump','uVA_B1','params_B1', 'bVA_B1', 'fRHS_B1', 'error_B1',...
%     'uVA_B2','params_B2', 'bVA_B2', 'fRHS_B2', 'error_B2',...
%     'uVA_B3','params_B3', 'bVA_B3', 'fRHS_B3', 'error_B3',...
%     'uVA_B4','params_B4', 'bVA_B4', 'fRHS_B4', 'error_B4');
