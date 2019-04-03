function rhs = GPE_nsoli_Stationary_1D(u,Nx,onehalfoodx2,V,V1, V2, dx, x)
global delta h sigma_X Ein tau0
U = reshape(u(1:Nx) + 1i*u([1:Nx]+Nx), Nx, 1);
mu = u(2*Nx+1);
%mu = 0; %real(u(2*Nx+1)) %abs(real(u(2*Nx+1))); %sqrt(u(2*Nx+1).*conj(u(2*Nx+1)))
MODU = U.*conj(U);
%dx = 0.1;

%%Fourth order
 ip=[1,1:Nx-1];
 im=[2:Nx,Nx];
 ipp=[1,1,1:Nx-2];
 imm=[3:Nx,Nx,Nx];

 %Second Order u''
 Um = [ U(Nx); U(1:Nx-1) ]; % periodic BCs
 Up = [ U(2:Nx); U(1) ];    % periodic BCs
 DU2 = Up-2*U+Um;
 
 %second Order phi''
 Vm = [ V(Nx); V(1:Nx-1) ]; % periodic BCs
 Vp = [ V(2:Nx); V(1) ];    % periodic BCs
 VU2 = Vp-2*V+Vm; 
 
 %first order u'
 DU = Up-Um;
 
 %first order v'
 VU = Vp-Vm;
 oodx = 0.5/dx;

% RHS = -onehalfoodx2*DU2 + (-mu + delta - MODU - 1i*(1+0.0)).*U + 1i*(V);  %test with adding a -i(1-alpha) term
 RHS = -onehalfoodx2*DU2 + (mu - MODU - 1i*(1+V2) +...
     V1.^2).*U - 2*1i*V1*oodx.*DU + 1i*(Ein);  
 rhs = [reshape(real(RHS),Nx,1); reshape(imag(RHS),Nx,1); mu];
 Uc = conj(U);
 Umc = [ Uc(Nx); Uc(1:Nx-1) ]; % periodic BCs
 Upc = [ Uc(2:Nx); Uc(1) ];    % periodic BCs
 DUc = Upc-Umc;
 rhs(2*Nx+1) = sum(-(1+V2).*MODU - V1.*oodx.*(DUc.*U + DU.*Uc) + 0.5*(U+conj(U)).*Ein)*dx; %0.5*(U+conj(U)).*Ein)*dx;% + 0.5*(U+conj(U)).*Ein)*dx;
% rhs(2*Nx+1) = sum(-(1+V2).*MODU - V1.*oodx.*(DUc.*U + DU.*Uc) + 0.5*(U+conj(U)).*Ein)*dx; %0.5*(U+conj(U)).*Ein)*dx;% + 0.5*(U+conj(U)).*Ein)*dx;
 %rhs(2*Nx+1) = sum(-(1+onehalfoodx2*VU2).*MODU - oodx.*VU.*oodx.*MODUDU +real(U) 0.5*(U+conj(U)).*Ein )*dx;% + 0.5*(U+conj(U)).*Ein)*dx;
 %rhs(2*Nx+1) = sum(0.5*(U+conj(U)).*Ein - MODU)*dx;
end