close all; clear all;
global allt allx 
global delta h sigma_X Ein tau0 uS

configLLERegular; 
load('RegularNCVA.mat')

cmax = [0.1, 0.5:0.5:20]';
betatau = [0.1, 0.5:0.5:20]'; 
NCVASurfPlot = ['NCVARegularMassPlot.mat'];
tt = [t0:dt:Tend];

NPotential = 280/2;
acc = 0.05;
NL = 1;
NR = 1601; 
tspan = [0 Tend];
tauStar = 2.5;

for i = 1:length(cmax);
    ct = cmax(i);
    NCVASaveName = ['NCVARegular_xf',num2str(ct),'.mat'];
    for k = 1:length(betatau);
        betat = betatau(k);
        [VAt,xVA]=ode45(@(t,x) NCVAVAF(t,x,betat,ct),tspan, p0);  
        stopsave = fix((length(VAt))/(saveframes-2));
        VAN(k) = length(VAt);
        for j = 1:length(VAt);
            if (j == 1)
                taut = 0;
                savenumber = 1; 
                uPt(k,savenumber) = 0;

                V = E(xi0, taut, sigma_X, h);
                V1 = dE(xi0, taut, sigma_X, h);
                V2 = ddE(xi0, taut, sigma_X, h);
                saveuVA(k,savenumber,:) = uVA;
                paramsNCVA(k,savenumber,:) = p0;
                Potential(k,savenumber,:) = V1.^2;

                uMass = uVA.*conj(uVA);
                result = round(taut/acc)*acc;
                indx = find(xi0 == result);
                Mass_In(savenumber) = sum(uMass(indx-NPotential:indx+NPotential))*dx;
                Mass_Out(savenumber) = sum(uMass(1:indx-NPotential-1))*dx + sum(uMass(indx+NPotential+1:end))*dx;
                MassIn(k,savenumber) =  Mass_In;
                MassOut(k,savenumber) = Mass_Out;
                MassTot = sum(uMass)*dx;
%                 plot(xi0, uVA.*conj(uVA), xi0, V1.^2)
%                 drawnow;
            end
            if (fix(j/stopsave)==j/stopsave)
                savenumber = savenumber +1; 
                t = VAt(j);
                uPt(k, savenumber) = t;
                if (t <= 2*tauStar)
                    taut = tau0 + (ct/2)*(tanh(betat*(t-tauStar))./tanh(betat*tauStar) + 1);
                else
                    taut = ct;
                end
                V1 = dE(xi0, taut, sigma_X, h);
                solp = xVA(j,:);
                uNCVA = (solp(1)).*exp(-0.5*(((xi0-solp(6))/solp(5)).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
                paramsNCVA(k,savenumber,:) = xVA(j,:);
                Potential(k,savenumber,:) = V1.^2;
                result = round(taut/acc)*acc;
                indx = find(floor(10000000000*(result-xi0))==0);
                uMass = uNCVA.*conj(uNCVA);
                Mass_In = sum(uMass(indx-NPotential:indx+NPotential))*dx;
                Mass_Out = sum(uMass(1:indx-NPotential-1))*dx + sum(uMass(indx+NPotential+1:end))*dx;
                MassIn(k,savenumber) =  Mass_In;
                MassOut(k,savenumber) = Mass_Out;
                saveuVA(k,savenumber,:) = uNCVA;
%                  plot(xi0, uNCVA.*conj(uNCVA), xi0, V1.^2)
%                  drawnow;pause(0.1);
            end  
            if (j == length(VAt))
                taut = ct;
                savenumber = savenumber + 1; 
                uPt(k,savenumber) = VAt(j);

                V1 = dE(xi0, taut, sigma_X, h);
                solp = xVA(j,:);
                uNCVA = (solp(1)).*exp(-0.5*(((xi0-solp(6))/solp(5)).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
                paramsNCVA(k,savenumber,:) = xVA(j,:);
                Potential(k,savenumber,:) = V1.^2;
                result = round(taut/acc)*acc;
                indx = find(floor(10000000000*(result-xi0))==0);
                uMass = uNCVA.*conj(uNCVA);
                Mass_In = sum(uMass(indx-NPotential:indx+NPotential))*dx;
                Mass_Out = sum(uMass(1:indx-NPotential-1))*dx + sum(uMass(indx+NPotential+1:end))*dx;
                MassIn(k,savenumber) =  Mass_In;
                MassOut(k,savenumber) = Mass_Out;
                saveuVA(k,savenumber,:) = uNCVA;
%                  plot(xi0, uNCVA.*conj(uNCVA), xi0, V1.^2)
%                   drawnow;pause(0.1);
                 saveN(k) = savenumber;
            end
        end
        MassTotal(k,:) = MassTot;
        QOut(k, i) = (MassOut(k,1) - MassOut(k,saveN(k)))/MassTot;
        QIn(k, i) = (MassIn(k,1) - MassIn(k,saveN(k)))/MassTot;    
    end
    save(NCVASaveName, 'saveuVA', 'MassIn','MassOut','Potential','uPt','ct','betatau','paramsNCVA','MassTotal','saveN','VAN')
end
save(NCVASurfPlot,'QOut','QIn','betatau','cmax')

