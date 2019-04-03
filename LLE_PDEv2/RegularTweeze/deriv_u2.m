%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dur = deriv_u2(u,Nx,ohdxx,eta,V,dx,V1,V2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global delta Ein
 um = [ u(Nx); u(1:Nx-1) ]; % periodic BCs
 up = [ u(2:Nx); u(1) ];    % periodic BCs
 DU2 = up-2*u+um;
 
 %second Order phi''
 vm = [ V(Nx); V(1:Nx-1) ]; % periodic BCs
 vp = [ V(2:Nx); V(1) ];    % periodic BCs
 DV2 = vp-2*V+vm; 
 
 %first order u'
 DU = up-um;
 
 %first order v'
 VU = vp-vm;
 oodx = 0.5/dx;
 
 oodx = 0.5/dx;
 ip=[1,1:Nx-1];
 im=[2:Nx,Nx];
 ipp=[1,1,1:Nx-2];
 imm=[3:Nx,Nx,Nx];
 MODU = u.*conj(u);
 
 %ux = oodx*(up-um);
 %dur = -1i.*(-ohdxx.*DU + (-1i*(1+0.0) - MODU + delta).*u + 1i.*V);

 dur = -1i.*(-ohdxx*DU2 + ( delta - MODU - 1i*(1+V2) +...
     (V1).^2).*u - 2*1i.*V1*oodx.*DU + 1i*(Ein));  %test with adding a -i(1-alpha) term


 return;

