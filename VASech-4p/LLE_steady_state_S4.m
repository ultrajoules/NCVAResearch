function RHS = LLE_steady_state_S4(x)
global delta h sigma_X Ein tau0
 warning('OFF','MATLAB:integral:MinStepSize');
 tt = [-100:0.1:100];
 FIa = @(t,x) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* (cosh((-t + x(4)) .* x(1)) + x(1) .* sinh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)).^ 2;
 Ia = integral(@(t)FIa(t,x),-50,50,'AbsTol',1e-15,'RelTol',1e-13);
 %integral %quadgT0 %trapz
 %FIa = @(t) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* (cosh((-t + x(4)) .* x(1)) + x(1) .* sinh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)).^ 2;
 %Ia = quadgk(FIa,-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 figure(1);
 plot(tt, FIa(tt,x))
 
 FId = @(t,x) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 Id = integral(@(t)FId(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %FId = @(t) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 %Id = quadgk(FId,-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIc = @(t,x) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* (-t + x(4)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 Ic = integral(@(t)FIc(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %FIc = @(t) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* (-t + x(4)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 %Ic = quadgk(FIc,-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIxi = @(t,x) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* (x(1) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* sinh((-t + x(4)) .* x(1)) + cos(x(2) + x(3) .* t - x(3) .* x(4)) .* cosh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)) .^ 2;
 Ixi = integral(@(t)FIxi(t,x),-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 %FIxi = @(t) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* (x(1) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* sinh((-t + x(4)) .* x(1)) + cos(x(2) + x(3) .* t - x(3) .* x(4)) .* cosh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)) .^ 2;
 %Ixi = quadgk(FIxi,-50,50,'AbsTol',1e-15,'RelTol',1e-13); 
 figure(2);
 plot(tt, FIxi(tt,x))
 warning('ON','MATLAB:integral:MinStepSize');

 % Q = (-iu+ delta*u + iS
 FA = -(2 * x(1)) + Id / 0.2e1;
 FD = ((2 * delta * x(1) + 2 * x(3) ^ 2 * x(1) + x(3) * Ic - Ia * x(1)) / x(1)) / 0.2e1;
 FC = -(x(3) * Id + Ixi) / x(1) / 0.2e1;
 FXi = ((4 * x(1) * x(3) + Ic) / x(1)) / 0.2e1;
 
 RHS =  [FA; FD; FC; FXi]; 
 
end