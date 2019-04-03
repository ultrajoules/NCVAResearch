function RHS = LLE_steady_state_S4(x)
global delta h sigma_X Ein tau0

 warning('ON','MATLAB:integral:MinStepSize');
 dt = 0.0001;
 tt = (-100:dt:100);
 x
 
  FIa = @(t,x) 0.2e1 .*  (Ein * exp(1i*(h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)))) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* (cosh((-t + x(4)) .* x(1)) + x(1) .* sinh((-t + x(4)) .* x(1))) .* sech((-t + x(4)) .* x(1)).^ 2; % ./ cosh((-t + x(4)) .* x(1)).^ 2;
%  Ia = integral(@(t)g(FIa(t,x)),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
  Ia = integral(@(t)FIa(t,x),-50,50,'AbsTol',1e-15,'RelTol',1e-13);
  %Ia = 2*delta-2*x(3)^2; 
  %Ia = sum(g(FIa(tt,x)))*dt;
 %integral %quadgT0 %trapz quadgk
 %FIa = @(t) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* (cosh((-t + x(4)) .* x(1)) + x(1) .* sinh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)).^ 2;
 %Ia = quadgk(FIa,-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 %figure(1);
 %plot(tt, FIa(tt,x))
 
 FId = @(t,x) 0.2e1 .*  (Ein * exp(1i*(h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)))) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) .*sech((-t + x(4)) .* x(1)); % ./ cosh((-t + x(4)) .* x(1));
 Id = integral(@(t)FId(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %FId = @(t) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 %Id = quadgk(FId,-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %Id = 4*x(1);
 % (Ein *  h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2))
 FIc = @(t,x) -0.2e1 .*  (Ein * exp(1i*(h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)))) .* (-t + x(4)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) .* sech((-t + x(4)) .* x(1)); %./ cosh((-t + x(4)) .* x(1));
 Ic = integral(@(t)FIc(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %FIc = @(t) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* (-t + x(4)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 %Ic = quadgk(FIc,-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %Ic = -4*x(1)*x(3);
 
 FIxi = @(t,x) -0.2e1 .* (Ein * exp(1i*(h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)))) .* x(1) .* (x(1) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* sinh((-t + x(4)) .* x(1)) + cos(x(2) + x(3) .* t - x(3) .* x(4)) .* cosh((-t + x(4)) .* x(1))) .* sech((-t + x(4)) .* x(1)) .^ 2; %./ cosh((-t + x(4)) .* x(1)) .^ 2;
 %Ixi = quadgk(@(t)g(FIxi(t,x)),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13); 
 Ixi = integral(@(t)FIxi(t,x),-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 %FIxi = @(t) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* (x(1) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* sinh((-t + x(4)) .* x(1)) + cos(x(2) + x(3) .* t - x(3) .* x(4)) .* cosh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)) .^ 2;
 %Ixi = quadgk(FIxi,-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 %figure(2);
 %plot(tt, FIxi(tt,x))
 %Ixi = sum(g(FIxi(tt,x)))*dt;
 %Ixi = -4*x(1)*x(3);
 warning('ON','MATLAB:integral:MinStepSize');
 Ixi
 Ia
 
 if isnan(Ixi)==1
     Ixi = 0;
 end 
 if isnan(Ia)==1
     Ia = 0;
 end 

 % Q = (-iu+ delta*u + iS
 FA = -(2 * x(1)) + Id / 0.2e1;
 FD = ((2 * delta * x(1) + 2 * x(3) ^ 2 * x(1) + x(3) * Ic - Ia * x(1)) / x(1)) / 0.2e1;
 FC = -(x(3) * Id + Ixi) / x(1) / 0.2e1;
 FXi = ((4 * x(1) * x(3) + Ic) / x(1)) / 0.2e1;
 
 RHS =  [FA; FD; FC; FXi]; 
 
 
end

  function y = g(z)
      y = z;
      y(isnan(y)) = 0;
  end