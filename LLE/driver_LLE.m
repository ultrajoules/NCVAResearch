%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %   
%  Julia Rossi                                        2/2016    %
%  LLE PDE solutions in time 
%                                                                %
%  driver_LLE.m                                                      %
%  Integrate the NLS for the solutions   %
%  using second-order central differencing in space and fourth-  %
%  order Runge Kutta in time. (note in paper t=z and x = T       %
%                                                                %  
% The NLS is given by:
%            uz - [-1 + 1i( |u|^2 - delta) - i*eta*uxx]u  -  E0(x) = 0        %
%                                                                %
% which is discretized into 
%       ut(t = n*k) = F(u) = [-1 + i( |u|^2 - delta) - i*eta*uxx]*u + E0(x)   
% where u0(x) = u0 + alpha*exp(-(x-tau0/beta)^2) 
%
%                                                                %
% There is a criteria on the discretization of space and time:
%      deltaT <= (less than or equal to) 
%                   2*sqrt(2)/(2+D) dx^2/onehalf
%    where D is the dimension of the problem
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

global allt allx 
global delta h sigma_X Ein tau0 uS
plotComparison = 0;
drawimg = 0;
PDEfit = 0;
Damping = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%File Saving names
NSOLIName = ['NSOLIResults/sol_',num2str(Ein),'.mat'];

config;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to add damping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Damping == 1)
    L = xmax;
    cutoff=L-10;
    wcut = (L-cutoff)/2;  % half-width of the cutoff region.
    xDamp=linspace(-L,L,N);

    mincut = 0.99;
    cutslope = 0.5; %0.5;

    damping_stepM=mincut+(1-mincut)*((1+tanh(cutslope*( xDamp+L-wcut)))/2);
    damping_stepP=mincut+(1-mincut)*((1+tanh(cutslope*(-xDamp+L-wcut)))/2);
    damping_step=damping_stepM.*damping_stepP;
    
    figure(4);clf
    plot(xDamp,damping_step)
    axis([-L L 0.9 1+(1-mincut)])
else 
    damping_step = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ii = 1:length(EinVec);

     Ein = EinVec(ii);
%     mu = delta;
%    for jj = 1:length(deltaVec);
                 
 %       delta = deltaVec(jj);
        for kk = 1:length(sigma_XVec);
             sigma_X = sigma_XVec(kk);
            for ll = 1:length(hVec);     
                h = hVec(ll);

                %mu= 1; %delta;
%                 if (h > 1)
%                    uini = uinit(xi0, 0, h,1, 1,0,tau0); %.*exp(1i*E(xi0,tau0));
%                 end

                Ein
                h
                sigma_X
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find initial conditions from NSOLI for PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = E(xi0, tau0, sigma_X, h); %Ein; %+ 0*E(xi0, tau0);
V1 = dE(xi0, tau0, sigma_X, h);
V2 = ddE(xi0, tau0, sigma_X, h);
%uSteadyState = -(1/6)*(-((108*1i)*Ein+12*sqrt(-12*1i+36*delta +(36*1i)*delta^2-12*delta^3-81*Ein^2))^(2/3)+12*1i-12*delta)/((108*1i)*Ein+12*sqrt(-12*1i+36*delta+(36*1i)*delta^2-12*delta^3-81*Ein^2))^(1/3)
%uSteadyState = +sqrt(delta) + 1i*Ein; 
%uSteadyState = Ein./(1+1i*(delta - uini.*conj(uini))); 0.3929 - 1i*0.7947; %sqrt(delta)+1i*Ein;
%uS = fsolve(@(uS) mySteadyState(uS,Ein, delta))
%x = fsolve(@(uS -uS+Ein./(1+1i*(delta - uS.*conj(uS))))
%uS = fsolve(Ein./(1+1i*(delta - uS.*conj(uS)))-uS); 0.3929 - 1i*0.7947; %sqrt(delta)+1i*Ein;
%uSteadyStateSol1 = -(-((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (0.2e1 / 0.3e1) - (12 * delta) + (12*1i)) * ((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (-0.1e1 / 0.3e1) / 6;
%uSteadyStateSol2 = (-((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (0.2e1 / 0.3e1) - (12 * delta) + (12*1i) + (1i) * sqrt(0.3e1) * ((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (0.2e1 / 0.3e1) + (-12*1i) * sqrt(0.3e1) * delta - (0.12e2 * sqrt(0.3e1))) * ((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (-0.1e1 / 0.3e1) / 12;
%uSteadyStateSol3 =  -(((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (0.2e1 / 0.3e1) + (12 * delta) + (-12*1i) + (1i) * sqrt(0.3e1) * ((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (0.2e1 / 0.3e1) + (-12*1i) * sqrt(0.3e1) * delta - (0.12e2 * sqrt(0.3e1))) * ((-108*1i) * Ein + 12 * sqrt(-(12 * delta ^ 3) + (36*1i) * (delta ^ 2) + (36 * delta) + (-12*1i) - 81 * Ein ^ 2)) ^ (-0.1e1 / 0.3e1) / 12;
%uini = (Ein/(sqrt(delta) +1) + 1i*Ein)*ones(N,1) + uini; %(-sqrt(delta)+1i*(Ein))*ones(N,1).*uini; %.*exp(1i*E(xi0,tau0));
%uini = uSteadyState+uini; %(-sqrt(delta)+1i*(Ein))*ones(N,1).*uini; %.*exp(1i*E(xi0,tau0));
% figure(10);
% plot(xi0, uini.*conj(uini))
%plot(xi0, uSteadyState.*uSteadyState, xi0, uSteadyStateSol1.*uSteadyStateSol1, xi0, uSteadyStateSol2.*uSteadyStateSol2, xi0, uSteadyStateSol3.*uSteadyStateSol3)
% ct = 0;
% betat = 0.01;
%  [allU, uPt] = cNLS1D(stopsave, saveframes, maxTiterations,N,dt,...
%             xi0,uini,ohdxx,eta,V,oos,allx,allt,damping_step,dx,E,dE,ddE,ct,betat);
% %pause;
% %upde = allU(length(allt),:);
% uini = allU(length(allt),:);

%NSOLI 
%mu = sigma_X-1; %delta;
mu = 3;
%uini = uinit(xi0, 0, 3, 1.5, 0.05,0,0.8, tau0); %.*exp(1i*E(xi0,tau0));
uini = uinit(xi0, 0, 3, 1.5, 0.05,0,0.8, tau0); %.*exp(1i*E(xi0,tau0));
uini = uini + (0.3929 - 0.7947i);
%uini = uini + (0.4 - 1i);

plot(xi0, uini.*conj(uini),xi0, real(uini), xi0, imag(uini),xi0, V1.^2);
mu =3;
%uini = u;
u0 = [reshape(real(uini),N,1); reshape(imag(uini),N,1); mu];
[uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,V1, V2,dx,xi0),1e-15*[1,1],sol_parms);

RHS = GPE_nsoli_Stationary_1D(uNSOLI,N,onehalfoodx2,V,V1, V2,dx,xi0);
figure(500); 
plot(1:2*N+1, RHS);
drawnow;

ubefore = reshape(uNSOLI(1:N) + 1i*uNSOLI([1:N]+N),N,1);
muNSOLI = uNSOLI(2*N+1);
mu = muNSOLI
theta = 0; %atan2(imag(ubefore), real(ubefore));
u = ubefore.*exp(-1i*theta);

figure(300);
%plot(xi0, u.*conj(u),xi0, real(u), xi0, imag(u), xi0, V.*conj(V), xi0,allU(length(allt),:).*conj(allU(length(allt),:)),  xi0, real(allU(length(allt),:)));
plot(xi0, u.*conj(u),xi0, real(u), xi0, imag(u), xi0, V, xi0, V1.^2); %.*conj(V));
title('Initial Condition after NSOLI');
%legend('|u|^2','Re(u)','Im(u)','|S(x)|^2','PDE |u|^2', 'PDE Re(u)');
legend('|u|^2','Re(u)','Im(u)','\phi(\tau)','\phi(\tau)^2');
%xlim([-20 20])
drawnow;

%save(NSOLIName, 'muNSOLI', 'allx', 'allt', 'u','V','ubefore');
%bifU = u(251).*conj(u(251));
uini= u;

testRHS = sum(RHS(1:2*N))/(2*N+1);
if (abs(RHS) < 1e-8)
%    stability1D;
    success(kk,ll) = 1;
     uS = fsolve(@(uS) -uS+Ein./(1+1i*(mu - uS.*conj(uS))), [1 1 1],optimoptions('fsolve','Display','iter', 'TolX', 1e-15,...
     'MaxIter', 100000, 'TolFun',1e-15));
     uSteadyState = ones(N,1)*uS(1);
 %    pks = findpeaks(uini.*conj(uini),'MinPeakHeight', uS(1).*conj(uS(1))+0.001);
 %    savepeaks(kk,ll,:) = pks;
%    eeSave(ii,jj,kk,ll,:) = ee;
%    vvSave(ii,jj,kk,ll,:,:) = vv;
else 
    mu = delta;
    uini = uinit(xi0, 0, h*3,sigma_X, 1,0,tau0); %.*exp(1i*E(xi0,tau0));
    u0 = [reshape(real(uini),N,1); reshape(imag(uini),N,1); mu];
    [uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,V1, V2,dx,xi0),1e-15*[1,1],sol_parms);

    RHS = GPE_nsoli_Stationary_1D(uNSOLI,N,onehalfoodx2,V,V1, V2,dx,xi0);
    figure(500); 
    plot(1:2*N+1, RHS);
    drawnow;
    ubefore = reshape(uNSOLI(1:N) + 1i*uNSOLI([1:N]+N),N,1);
    muNSOLI = uNSOLI(2*N+1);
    mu = muNSOLI
    testRHS = sum(RHS(1:2*N))/(2*N+1);
    if (abs(RHS) < 1e-8)
        success(kk,ll) = 1;
        uS = fsolve(@(uS) -uS+Ein./(1+1i*(mu - uS.*conj(uS))), [1 1 1],optimoptions('fsolve','Display','iter', 'TolX', 1e-15,...
          'MaxIter', 100000, 'TolFun',1e-15));
        uSteadyState = ones(N,1)*uS(1);
  %      pks = findpeaks(uini.*conj(uini),'MinPeakHeight', uS(1).*conj(uS(1))+0.001);
  %      savepeaks(kk,ll,:) = pks;
    else
        success(kk,ll)=0;
    end
end


uS
mu
 Vm = [ V(N); V(1:N-1) ]; % periodic BCs
 Vp = [ V(2:N); V(1) ];    % periodic BCs
uSave(kk,ll,:) = u;
%uSaveSteadyState(ii,kk,ll,:) = uS;
nsoliError(kk,ll) = ierr;
phiSave(kk,ll,:) = V;
phidotSave(kk,ll,:) = V1;
muDelta(kk,ll) = muNSOLI;
            end
        end
end
    
uS = uS(1);
delta = mu;

save('FatLLE','uini','u','uSteadyState','uS','mu','delta','N','xi0','Ein','sigma_X','h')

%Now need to find the steady state solution without no soliton 

pause; 
uiniNOSol = uSteadyState;

u0NoSol = [reshape(real(uiniNOSol),N,1); reshape(imag(uiniNOSol),N,1); mu];
[uNoSol,it_hist,ierr,u_hist]=nsoli(u0NoSol,@(uNoSol_)GPE_nsoli_Stationary_1DSS(uNoSol_,N,onehalfoodx2,V,V1, V2,dx,xi0),1e-15*[1,1],sol_parms);

RHS = GPE_nsoli_Stationary_1D(uNoSol,N,onehalfoodx2,V,V1, V2,dx,xi0);
figure(500); 
plot(1:2*N+1, RHS);
drawnow;

uNS = reshape(uNoSol(1:N) + 1i*uNoSol([1:N]+N),N,1);
muNS = uNoSol(2*N+1);

figure(300);
plot(xi0, uNS.*conj(uNS),xi0, real(uNS), xi0, imag(uNS), xi0, V, xi0, V1.^2); %.*conj(V));
title('Initial Condition after NSOLI');
legend('|u|^2','Re(u)','Im(u)','\phi(\tau)','\phi(\tau)^2');
drawnow;


 %save('SkinnyLLE','uini','u','uSteadyState','uS','mu','delta','N','xi0','Ein','sigma_X','h')
 %save('RegularLLE','uini','u','uSteadyState','uS','mu','delta','N','xi0','Ein','sigma_X','h')
 save('FatLLENoSol','uNS','muNS')
% 
% for k =2;
% for i = 1:length(uPt);
% test(:) = saveU(k,i,:);
% V(:) = Potential(k,i,:);
% plot(xi0, test.*conj(test), xi0, V)
% drawnow;
% end
% end
 
 
% uini= u;
% %pause;
% delta = muNSOLI;
% cmax = [0.5:0.5:5];
% betatau = [0.05,0.1:0.1:1,1.5,2,2.5,3]; %[0.1:0.1:1];
% uSS = uSteadyState.*conj(uSteadyState);
% for i = 1:length(cmax);
%     for k = 1:length(betatau);
%         ct = cmahx(i);
%         betat = betatau(k);
%         [allU, uPt, Mass_In, Mass_Out] = cNLS1D(stopsave, saveframes, maxTiterations,N,dt,...
%             xi0,uini,ohdxx,eta,V,oos,allx,allt,damping_step,dx,E,dE,ddE,ct,betat,uSS);
%         QOut(k, i) = (Mass_Out(1) - Mass_Out(end))/Mass_Out(1)
%         QIn(k, i) = (Mass_In(1) - Mass_In(end))/Mass_In(1)
%     end
% end
% tSave = uPt;
% 
% %This should be done in the first loop!!!! 
% for i = 1:length(cmax);
%     for k = 1:length(betatau);
%         if (k ==1)
%             Utest(i,k,:) = USave_b005(i,end,:);
%         elseif (k==2)
%             Utest(i,k,:) = USave_b01(i,end,:);
%         elseif (k==3)
%             Utest(i,k,:) = USave_b02(i,end,:);
%         elseif k==4 
%             Utest(i,k,:) = USave_b03(i,end,:);
%         elseif (k==5)
%             Utest(i,k,:) = USave_b04(i,end,:);
%         elseif (k==6)
%             Utest(i,k,:) = USave_b05(i,end,:);
%         elseif (k==7)
%            Utest(i,k,:) =  USave_b06(i,end,:);
%         elseif (k==8)
%             Utest(i,k,:) = USave_b07(i,end,:);
%         elseif (k==9)
%             Utest(i,k,:) = USave_b08(i,end,:);
%         elseif (k==10)
%             Utest(i,k,:) = USave_b09(i,end,:);
%         elseif (k==11)
%             Utest(i,k,:) = USave_b10(i,end,:);
%         elseif (k==12)
%             Utest(i,k,:) = USave_b15(i,end,:);
%         elseif (k==13)
%             Utest(i,k,:) = USave_b20(i,end,:);
%         elseif (k==14)
%             Utest(i,k,:) = USave_b25(i,end,:);
%         elseif (k==15)
%             Utest(i,k,:) = USave_b30(i,end,:); 
%         end
%     end
% end
%         
% 
% ep = 2;
% Uo = max(u.*conj(u));
% j = 0;
% l = 0; 
% for i = 1:length(cmax); 
%     for k = 1:length(betatau);
%         Ut = squeeze(Utest(i,k,:)); 
%         max(Ut.*conj(Ut))
%         if (max(Ut.*conj(Ut))) > Uo-ep  & (max(Ut.*conj(Ut))) < Uo+ep 
%             statusTweeze(i,k) = 1; 
%             j = j+1;
%             TweezeP(j,:) = [cmax(i), betatau(k)];
%         else
%             statusNoTweeze(i,k) = 1; 
%             l = l+1;
%             NoTweezeP(l,:) = [cmax(i), betatau(k)];
%         end
%     end
% end
% 
% [~,idb] = max(statusNoTweeze(:,sum(statusNoTweeze)>0));
% idc = [4:1:15];
% plot(betatau(idc), cmax(idb),'LineWidth',4)
% xlabel('$\beta$','Interpreter', 'LaTex','Fontsize',18)
% ylabel('$c_\infty$','Interpreter', 'LaTex','Fontsize',18)
% 
% %For plotting - extract the lowest value of c for each beta
% % NoTweezeSort = sortrows(NoTweezeP,2)
% % plotLine = [1,3,7,13,20,27,35,43,51,60,69,79];
% % for i = 1:length(plotLine)
% %     plotNoTweezeLine(i,:) = NoTweezeSort(plotLine(i),:);
% % end
%         
% 
%     
% 
% 
% 
% % jj = 0;
% % pp = 0;
% % for ii = 9:length(EinVec);
% %     
% %      Ein = EinVec(ii);
% % %     mu = delta;
% % %    for jj = 1:length(deltaVec);
% %                  
% %  %       delta = deltaVec(jj);
% %         for kk = 1:length(sigma_XVec);
% %              sigma_X = sigma_XVec(kk);
% %             for ll = 1:length(hVec);     
% %                 h = hVec(ll);
% %                 Ein
% %                 h
% %                 sigma_X
% %                 if (success(ii,kk,ll) == 0 )
% %                     if (ll >= 2)
% %                         uini(:) = uSave(ii,kk,ll-1,:);
% %                         mu = muDelta(ii,kk,ll-1);
% %                         uini = uini';
% %                     else 
% %                         uini = uinit(xi0, 0, h*2,sigma_X, h,0,tau0); 
% %                         mu = muDelta(ii,kk,ll+1); 
% %                     end
% %                     V = E(xi0, tau0, sigma_X, h);
% %                     u0 = [reshape(real(uini),N,1); reshape(imag(uini),N,1); mu];
% %                     [uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,dx,xi0),1e-15*[1,1],sol_parms);
% % 
% %                     RHS = GPE_nsoli_Stationary_1D(uNSOLI,N,onehalfoodx2,V,dx,xi0);
% %                     figure(500); 
% %                     plot(1:2*N+1, RHS);
% %                     drawnow;
% % 
% %                     ubefore = reshape(uNSOLI(1:N) + 1i*uNSOLI([1:N]+N),N,1);
% %                     muNSOLI = uNSOLI(2*N+1);
% %                     mu = muNSOLI;
% %                     theta = 0; %atan2(imag(ubefore), real(ubefore));
% %                     u = ubefore.*exp(-1i*theta);
% % 
% %                     figure(300);
% %                     %plot(xi0, u.*conj(u),xi0, real(u), xi0, imag(u), xi0, V.*conj(V), xi0,allU(length(allt),:).*conj(allU(length(allt),:)),  xi0, real(allU(length(allt),:)));
% %                     plot(xi0, u.*conj(u),xi0, real(u), xi0, imag(u), xi0, V); %.*conj(V));
% %                     title('Initial Condition after NSOLI');
% %                     %legend('|u|^2','Re(u)','Im(u)','|S(x)|^2','PDE |u|^2', 'PDE Re(u)');
% %                     legend('|u|^2','Re(u)','Im(u)','\phi(\tau)');
% %                     xlim([-20 20])
% %                     drawnow;
% % 
% %                     testRHS = sum(RHS(1:2*N))/(2*N+1);
% %                     if (abs(RHS) < 1e-8)
% %                     %    stability1D;
% %                         success(ii,kk,ll) = 1;
% %                     %    eeSave(ii,jj,kk,ll,:) = ee;
% %                     %    vvSave(ii,jj,kk,ll,:,:) = vv;
% %                     else 
% %                         success(ii,kk,ll) = 0;
% %                     end
% %                     Vm = [ V(N); V(1:N-1) ]; % periodic BCs
% %                     Vp = [ V(2:N); V(1) ];    % periodic BCs
% %                     uSave(ii,kk,ll,:) = u;
% %                     %uSaveSteadyState(ii,kk,ll,:) = uS;
% %                     nsoliError(ii,kk,ll) = ierr;
% %                     phiSave(ii,kk,ll,:) = V;
% %                     phidotSave(ii,kk,ll,:) = oodx*(Vp-Vm);
% %                     muDelta(ii,kk,ll) = muNSOLI;
% %                     savePeaks(ii,kk,ll) = NaN;
% %                 end
% %                 
% %                 if (success(ii,kk,ll) == 1)
% %                     delta = muDelta(ii,kk,ll);
% %                     uS = fsolve(@(uS) -uS+Ein./(1+1i*(delta - uS.*conj(uS))), [1 1 1],optimoptions('fsolve','Display','iter', 'TolX', 1e-15,...
% %                         'MaxIter', 100000, 'TolFun',1e-15));
% %                     uSteadyState = ones(N,1)*uS(1);
% %                     uini(:) = uSave(ii,kk,ll,:);
% %                     pks = findpeaks(uini.*conj(uini),'MinPeakHeight', uS(1).*conj(uS(1))+0.001);
% %                     if (isempty(pks) == 1)
% %                         savePeaks(ii,kk,ll) = 0;
% %                     else 
% %                         if (length(pks) ~= 0 ) 
% %                             savePeaks(ii,kk,ll) = length(pks);
% %                             if (length(pks) == 1);
% %                                  jj = jj + 1;
% %                                  E1Peak(jj,:) = Ein;
% %                                  E1VSigma(jj,:) = [Ein, sigma_X];
% %                                  E1Vh(jj,:) = [Ein, h];
% %                                  Sigma1VhVdelta(jj,:) = [sigma_X, h, delta];
% %                                  E1VSigmaVhVdelta(jj,:) = [Ein, sigma_X, h, delta];
% %                                  u1Peak(jj,:) = uSave(ii,kk,ll,:);
% %                                  u1SteadyState(jj,:) = uS(1);
% %                             else
% %                                  pp = pp + 1;
% %                                  EMPeak(pp,:) = Ein;
% %                                  EMVSigma(pp,:) = [Ein, sigma_X];
% %                                  EMVh(pp,:) = [Ein, h];
% %                                  SigmaMVhVdelta(pp,:) = [sigma_X, h, delta];
% %                                  EMVSigmaVhVdelta(pp,:) = [Ein, sigma_X, h, delta];
% %                                  uMPeak(pp,:) = uSave(ii,kk,ll,:);  
% %                                  uMSteadyState(pp,:) = uS(1);
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% % end
% 
%            
%                 
%                     
%                     
%                     
%                     
% 
% 
% %save(SSBName, 'PP', 'bifU')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fit PDE to ansatz 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = u-uSteadyState;
% %PDEFitting_S4;
% PDEFitting_G6;
% %PDEFitting_G4Tau;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Find NCVA solution for the variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uS = uSS(1);
% %NCVA_LLE_S4;
% NCVA_LLE_G6;
% %NCVA_LLE_G4Tau;
% %plot(xi0, uVA.*conj(uVA), xi0, u.*conj(u)); %, xi0, uppde.*conj(uppde)) %/(u.^2))
%     
% if (drawimg == 1) 
%     figure(2);
%     for jj=1:length(allt);
%         plot(xi0,allU(jj,:).*conj(allU(jj,:)), xi0, uVA(jj,:).*conj(uVA(jj,:))); %'o', xi0, uVANC(jj,:).*conj(uVANC(jj,:)));
%         title(['T = ', num2str(allt(jj))]);
%         xlabel('$x$','Interpreter', 'LaTex','Fontsize',18)
%         ylabel('$|u(x,t)|^2$','Interpreter', 'LaTex','Fontsize',18)
%         legend('PDE','NCVA'); %'VA', 'NCVA');
%         axis([-xmax xmax 0 maxA+0.1]);
%         drawnow;
% %       pause(0.01);
%     end
% end
% 
% 
% if (plotComparison == 1)
% 
%     %Plot the parameters of ODE solution vs. PDE solution 
%     % over time to see how they compare!
%     %AA = allANC.*exp(-imag(allBNC));
%     figure(30);
%     subplot(3,2,[1,2])
%     plot(x, allU(1,:).*conj(allU(1,:)),'b-.',x,uVA(1,:).*conj(uVA(1,:)),'r.',...
%      x, allU(saveframes+1,:).*conj(allU(saveframes+1,:)),'g-o', x, uVA(saveframes+1,:).*conj(uVA(saveframes+1,:)),'k.',...
%      'LineWidth', 4, 'MarkerSize',10);
%     xlim([x0-10, -x0+10]);
%     %ylim([0,1]);
%     title(['\bf{PDE VS ODE Dynamics, $\chi =$ ', num2str(chi),', $\sigma = $', num2str(sigma),'}'],'Interpreter', 'LaTex','Fontsize',36);
%     xlabel('x','Interpreter', 'LaTex','Fontsize',24);
%     ylabel('$|u(x,t)|^2$','Interpreter', 'LaTex','Fontsize',24);
%     legend('PDE: t=0','NCVA: t=0', ['PDE: t=', num2str(allt(saveframes+1))], ['NCVA: t=', num2str(allt(saveframes+1))]);
%     h = legend;
%     set(h,'Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
% 
%     subplot(3,2,3)
%     plot(allt, real(pp(:,1)),'b-', allt, allA,'r.','Linewidth',4,'MarkerSize',10)
%     %ylim([0.8, 1.2]);
%     xlim([0 Tend]);
%     ylim([0.9 1.1]);
%     %ylim([0.9 2.1])
%     xlabel('t','Interpreter','LaTex','Fontsize',24);
%     ylabel('$a$','Interpreter','LaTex','Fontsize',24);
%     title('Parameter Comparison for $a(t)$','Interpreter','LaTex','Fontsize',24);
%     legend('PDE','NCVA');
%     h = legend;
%     set(h,'Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     %axis tight;
% 
%     subplot(3,2,4)
%     plot(allt, real(pp(:,2)),'b-', allt, allB,'r.','Linewidth', 4,'MarkerSize',10)
%     %plot(allt, real(pp(:,2)),'b-o', allt, allB,'ro',allt, 0.5*allt'.*(allANC.*exp(imag(allBNC)).^2 + (allCNC).^2),'k+')
%     xlim([0 Tend]);
%     xlabel('$t$','Interpreter','LaTex','Fontsize',24);
%     ylabel('$b$','Interpreter','LaTex','Fontsize',24);
%     title('Parameter Comparison for $b(t)$','Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     %legend('PDE','ODE');
%     axis tight;
% 
% 
%     subplot(3,2,5)
%     plot(allt, real(pp(:,3)),'b-', allt, allC,'r.','Linewidth',4,'MarkerSize',10)
%     %ylim([-0.1, 0.3]);
%     ylim([cc-0.2, cc+0.2]);
%     xlim([0 Tend]);
%     xlabel('$t$','Interpreter','LaTex','Fontsize',24);
%     ylabel('$c$','Interpreter','LaTex','Fontsize',24);
%     title('Parameter Comparison for $c(t)$','Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     %legend('PDE','ODE');
% 
% 
%     subplot(3,2,6)
%     plot(allt, real(pp(:,4)), 'b-',allt, allXi,'r.','Linewidth',4,'MarkerSize',10)
%     xlim([0 Tend]);
%     xlabel('$t$','Interpreter','LaTex','Fontsize',24);
%     ylabel('$\xi$', 'Interpreter','LaTex','Fontsize',24);
%     title('Parameter Comparison for $\xi(t)$', 'Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     %legend('PDE','ODE');
%     axis tight;
% 
% 
%     figure(31)
%     %subplot(3,1,1)
%     %title(['PDE VS ODE Dynamics, $\epsilon =$ ', num2str(ep)],'Interpreter', 'LaTex','Fontsize',18);
%     %subplot(3,1,1)
%     plot(allt, real(pp(:,1)),'b-', allt, allA,'r.','Linewidth',4,'MarkerSize',10)
%     xlabel('$t$','Interpreter','LaTex','Fontsize',24);
%     ylabel('$a$','Interpreter','LaTex','Fontsize',24);
%     ylim([0.9, 1.1]);
%     title({...
%         ['\bf{PDE VS ODE Dynamics, $\chi =$ ', num2str(chi),', $\sigma = $', num2str(sigma), '}']
%         'Parameter Comparison for $a(t)$'
%         },'Interpreter','LaTex','Fontsize',36);
%     legend('PDE','NCVA');
%     h = legend;
%     set(h,'Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     %axis tight;
% 
% 
%     %%%relative error of ODE to PDE for a 
%     relerrPVA = abs((real(pp(:,1)) - allA'))./real(pp(:,1));
%     %relerrNCVA = abs((real(pp(:,1)) - AA))./real(pp(:,1));
%     figure(35)
%     plot(allt, relerrPVA, 'bs','Linewidth',2)
%     title('\bf{Relative Error in ODE Parameter $a(t)$}','Interpreter','LaTex','Fontsize',36)
%     xlabel('$t$','Interpreter','LaTex','Fontsize',24);
%     ylabel('Relative Error','Interpreter','LaTex','Fontsize',24);
%     %legend('P-VA','NC-VA');
%     %h = legend;
%     set(h,'Interpreter','LaTex','Fontsize',24);
%     set(gca,'LineWidth',2,'Fontsize',24);
%     % 
% 
%    
% end

