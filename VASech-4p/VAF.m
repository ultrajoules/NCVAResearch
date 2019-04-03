function VAx = VAF(tt, x)
VAx = zeros(4,1);  %output must be a column vector

global delta h sigma_X Ein tau0
 warning('OFF','MATLAB:integral:MinStepSize');

 FIa = @(t,x) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* (cosh((-t + x(4)) .* x(1)) + x(1) .* sinh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)).^ 2;
 Ia = integral(@(t)FIa(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %integral %quadgT0 %trapz
 FId = @(t,x) 0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 Id = integral(@(t)FId(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIc = @(t,x) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* (-t + x(4)) .* x(1) .* cos(x(2) + x(3) .* t - x(3) .* x(4)) ./ cosh((-t + x(4)) .* x(1));
 Ic = integral(@(t)FIc(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIxi = @(t,x) -0.2e1 .* (Ein + h .* exp(-(t - tau0).^ 2 ./ sigma_X ^ 2)) .* x(1) .* (x(1) .* sin(x(2) + x(3) .* t - x(3) .* x(4)) .* sinh((-t + x(4)) .* x(1)) + cos(x(2) + x(3) .* t - x(3) .* x(4)) .* cosh((-t + x(4)) .* x(1))) ./ cosh((-t + x(4)) .* x(1)) .^ 2;
 Ixi = integral(@(t)FIxi(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13); 
 warning('ON','MATLAB:integral:MinStepSize');

 FA = -(2 * x(1)) + Id / 0.2e1;
 FD = ((2 * delta * x(1) + 2 * x(3) ^ 2 * x(1) + x(3) * Ic - Ia * x(1)) / x(1)) / 0.2e1;
 FC = -(x(3) * Id + Ixi) / x(1) / 0.2e1;
 FXi = ((4 * x(1) * x(3) + Ic) / x(1)) / 0.2e1;
 
 RHS =  [FA; FD; FC; FXi]; 
 VAx(1) = FA;
 VAx(2) = FD;
 VAx(3) = FC;
 VAx(4) = FXi;
 
end

