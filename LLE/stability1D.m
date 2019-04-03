%
% This code computes stability for GPE
%       i u_t = -(onehalf) (u_xx) + g*|u|^2 u + V(x) u +i(Chi(x) -
%       sigma*|u|^2)*u
%  g = alpha = 1 is defocusing
%  g = alpha =-1 is focusing
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   allAumf = [];
    neigs=100; % # of evals to compute. neigs=0 => FULL spectrum
    allee=[];
    
 fprintf('setting up matrices\n');
 LAP=onehalfoodx2*LapMatrix1d_0bc(N);
 oodx = 0.5/dx;
 
 uu=u.*u;
 u2=u.*conj(u);
 
  %second Order phi''
%  Vm = [ V(N); V(1:N-1) ]; % periodic BCs
%  Vp = [ V(2:N); V(1) ];    % periodic BCs
%  VU2 = onehalfoodx2*(Vp-2*V+Vm); 
%   %first order phi'
%  VU = oodx*(Vp-Vm);
  vv = V1.^2;
%  VU=+diag(2*1i*V(1:N-1),1)-diag(2*1i*Vm(2:N),-1);
%  VU= oodx*VU;
%  CD = oodx*VU;
 centralDiff = CentralDiff1d_0bc(N);
 %CD = 2*1i*V1*centralDiff;
 CD = +diag(2*1i*V1(1:N-1),1)-diag(2*1i*V1(2:N),-1);
 
 
 %M11=diag(-muNSOLI+V+2*(g-1i*sigma)*u2+1i*Chi);
 %M12=diag((g-1i*sigma)*uu);

 M11=diag(-2*g*u2-(delta-1i*1)-vv+1i*V2);
 M12=diag(-(g)*uu);

 J=(1i*[LAP+M11+CD,M12;-conj(M12),-conj(LAP+M11+CD)]);

 neigs = N;
 if(exist('neigs')&neigs>0)
  optionseigs.disp=0;
  [vvv,eee]=eigs(sparse(J),neigs,'lr',optionseigs);  %0.5 for anchor so alwyas check or put 'lr'-largest to the real
 else
  [vvv,eee]=eig(J);
 end

ee=diag(eee);                        %eigenvalues
vv=vvv(1:N,:)+conj(vvv(N+1:end,:));  %eigenvectors
[aa,bb]=(sort(real(ee)));

%plot_eigen

