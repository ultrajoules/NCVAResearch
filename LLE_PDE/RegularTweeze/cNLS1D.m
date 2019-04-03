function [allU, uPt, Mass_In, Mass_Out, potential] = cNLS1D(stopsave, saveframes, maxTiterations,N,dt,xi0,...
    uini,ohdxx,eta,oos,allx,allt,damping_step,dx,E,dE,ddE,ct,betat,uSS)
global delta Ein sigma_X h tau0 uS
allU = zeros(saveframes+1,N);
Usave = zeros(length(xi0), saveframes+1);
savex = allx';
savet = allt';

tauStar = 2.5;
NPotential = 152/2;
acc = 0.05;
NL = 1;
NR = 1601; 
%define indices as row vectors
%periosdic BCS
i0 = [1:N];
ip = [i0(2:N), i0(1)];
im = [i0(N), i0(1:N-1)];

%initialize%
%U = zeros(N,1);
U = uini;

%figure('Renderer','zbuffer')
 %figure(1);
 %plot(xi0(NL:NR), real(U(NL:NR)), 'g', xi0(NL:NR), imag(U(NL:NR)), 'r', xi0(NL:NR), U(NL:NR).*conj(U(NL:NR)), 'b'); 
 %xlim([-20 20])
 %drawnow;
 
 %figure(3);
 %plot(allt,  tau0 + (ct/2)*(tanh(betat*(allt-tauStar))./tanh(betat*tauStar) + 1));
 %drawnow;
 %figure(4);
 
t = 0;
taut = tau0;
savenumber = 0;
savenumber = savenumber + 1;
allU(savenumber,:) = U;
Usave(:,savenumber) = U';
uPt(1) = t;
V = E(xi0, taut, sigma_X, h);
V1 = dE(xi0, taut, sigma_X, h);
V2 = ddE(xi0, taut, sigma_X, h);
uMass = U.*conj(U);
result = round(taut/acc)*acc;
indx = find(xi0 == result);
Mass_In(savenumber) = sum(uMass(indx-NPotential:indx+NPotential)-uSS(indx-NPotential:indx+NPotential))*dx;
Mass_Out(savenumber) = sum(uMass(1:indx-NPotential-1)-uSS(1:indx-NPotential-1))*dx + sum(uMass(indx+NPotential+1:end) - uSS(indx+NPotential+1:end))*dx;
for j = 1:maxTiterations;
    t = t+dt;
    taut = tau0 + (ct/2)*(tanh(betat*(t-tauStar))./tanh(betat*tauStar) + 1); 
    V = E(xi0, taut, sigma_X, h);
    V1 = dE(xi0, taut, sigma_X, h);
    V2 = ddE(xi0, taut, sigma_X, h);
    uN = ode_rk_dissipation(U,N,ohdxx,eta,V,dt,oos,dx,V1,V2);
    U = uN; 
    if (fix(j/stopsave)==j/stopsave)
        savenumber = savenumber + 1;
        uPt(savenumber) = t;
        allU(savenumber,:) = U(:);
        Usave(:,savenumber) = U';
        uMass = U.*conj(U);
        %plot(xi0(NL:NR), real(U(NL:NR)), xi0(NL:NR), imag(U(NL:NR)), xi0(NL:NR), U(NL:NR).*conj(U(NL:NR)),xi0(NL:NR), V1(NL:NR).^2);
        %xlim([-20 20])
        %drawnow;
        potential(savenumber,:) = V1.^2;
        result = round(taut/acc)*acc;
        indx = find(floor(10000000000*(result-xi0))==0);
        Mass_In(savenumber) = sum(uMass(indx-NPotential:indx+NPotential) - uSS(indx-NPotential:indx+NPotential))*dx;
        Mass_Out(savenumber) = sum(uMass(1:indx-NPotential-1) - uSS(1:indx-NPotential-1))*dx + sum(uMass(indx+NPotential+1:end) - uSS(indx+NPotential+1:end))*dx;
    end  
end

end