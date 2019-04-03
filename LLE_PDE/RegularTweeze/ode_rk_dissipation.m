

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rkZn = ode_rk_dissipation(rku,Nx,ohdxx,eta,V,dt,oos,dx,V1,V2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global delta Ein 
 k1 = dt*deriv_u2(rku       ,Nx,ohdxx,eta,V,dx,V1,V2);
 k2 = dt*deriv_u2(rku+0.5*k1,Nx,ohdxx,eta,V,dx,V1,V2);
 k3 = dt*deriv_u2(rku+0.5*k2,Nx,ohdxx,eta,V,dx,V1,V2);
 k4 = dt*deriv_u2(rku+    k3,Nx,ohdxx,eta,V,dx,V1,V2);

 rkZn = rku + oos*(k1 + 2*k2 + 2*k3 + k4);

 return;
