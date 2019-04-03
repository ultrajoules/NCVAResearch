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

configLLEFat; 
load('FatLLE.mat');
load('FatLLENoSol.mat');

cmax = [0.1, 0.5:0.5:20]';
betatau = [0.1, 0.5:0.5:20]'; 
uNS2 = uNS(end)*ones(1600,1);
uNNS = [uNS2; uNS];
uSS = uNNS.*conj(uNNS);
PDESurfPlot = ['LLEFatMassPlot.mat'];

for i = 1:length(cmax);
    ct = cmax(i);
    PDESaveName = ['LLEFat_xf',num2str(ct),'.mat'];
    for k = 1:length(betatau);
        betat = betatau(k);
        [allU, uPt, Mass_In, Mass_Out, potential,MassTot] = cNLS1D(stopsave, saveframes, maxTiterations,N,dt,...
            xi0,uini,ohdxx,eta,oos,allx,allt,damping_step,dx,E,dE,ddE,ct,betat,uSS, xNS);
        QOut(k, i) = (Mass_Out(1) - Mass_Out(end))/MassTot;
        QIn(k, i) = (Mass_In(1) - Mass_In(end))/MassTot;  
        saveU(k,:,:) = allU;
        MassIn(k,:) =  Mass_In;
        MassOut(k,:) = Mass_Out;
        Potential(k,:,:) = potential;
        MassTotal(k,:) = MassTot;
    end
    save(PDESaveName, 'saveU', 'MassIn','MassOut','Potential','uPt','ct','betatau','MassTotal')
end
save(PDESurfPlot,'QOut','QIn','betatau','cmax')

%driver_FatMassPlots