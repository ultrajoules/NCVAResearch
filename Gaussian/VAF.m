function VAx = VAF(tt, x)
VAx = zeros(6,1);  %output must be a column vector

global delta T0 P
 warning('OFF','MATLAB:integral:MinStepSize');

 FIa = @(t,x) 0.2e1 * exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 * t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * sqrt(P) .* sin((x(4) * t.^2 - 2 * x(4) * t * x(6) + x(4) * x(6) ^ 2 + x(3) * t - x(3) * x(6) + x(2)));
 Ia = integral(@(t)FIa(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 %integral %quadgT0 %trapz
 FIb = @(t,x) 0.2e1 * sqrt(P) * exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 * t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * x(1) .* cos((x(4) * t.^2 - 2 * x(4) * t * x(6) + x(4) * x(6) ^ 2 + x(3) * t - x(3) * x(6) + x(2)));
 Ib = integral(@(t)FIb(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIc = @(t,x) -0.2e1 * sqrt(P) * exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 * t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * x(1) .* (-t + x(6)) .* cos((x(4) * t.^2 - 2 * x(4) * t * x(6) + x(4) * x(6) ^ 2 + x(3) * t - x(3) * x(6) + x(2)));
 Ic = integral(@(t)FIc(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FId = @(t,x) 0.2e1 * sqrt(P) * exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 * t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * x(1) .* ((-t + x(6)) .^ 2) .* cos((x(4) * t.^2 - 2 * x(4) * t * x(6) + x(4) * x(6) ^ 2 + x(3) * t - x(3) * x(6) + x(2)));
 Id = integral(@(t)FId(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIs = @(t,x) 0.2e1 * ((-t + x(6)) .^ 2) * x(1) .* exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 * t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * sqrt(P) .* sin((x(4) * t.^2 - 2 * x(4) * t * x(6) + x(4) * x(6) ^ 2 + x(3) * t - x(3) * x(6) + x(2))) / (x(5) ^ 3);
 Is = integral(@(t)FIs(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);

 FIxi = @(t,x) -0.2e1 * sqrt(P) .* exp(-((T0 ^ 2 * t.^2 - 2 * T0 ^ 2 .* t * x(6) + T0 ^ 2 * x(6) ^ 2 + 2 * t.^2 * x(5) ^ 2) / x(5) ^ 2 / T0 ^ 2) / 0.2e1) * x(1) .* (0.2e1 .* cos((x(4) * t.^2 - 2 * x(4) .* t * x(6) + x(4) * x(6) ^ 2 + x(3) .* t - x(3) * x(6) + x(2))) .* (x(5) ^ 2) * x(4) .* t - 0.2e1 .* cos((x(4) * t.^2 - 2 * x(4) .* t * x(6) + x(4) * x(6) ^ 2 + x(3) .* t - x(3) * x(6) + x(2))) .* (x(5) ^ 2) * x(4) * x(6) + cos((x(4) * t.^2 - 2 * x(4) .* t * x(6) + x(4) * x(6) ^ 2 + x(3) .* t - x(3) * x(6) + x(2))) .* (x(5) ^ 2) * x(3) - sin((x(4) * t.^2 - 2 * x(4) .* t * x(6) + x(4) * x(6) ^ 2 + x(3) .* t - x(3) * x(6) + x(2))) .* t + sin((x(4) * t.^2 - 2 * x(4) .* t * x(6) + x(4) * x(6) ^ 2 + x(3) .* t - x(3) * x(6) + x(2))) * x(6)) / (x(5) ^ 2);
 Ixi = integral(@(t)FIxi(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13); 
 warning('OFF','MATLAB:integral:MinStepSize');

 % Q = (-iu+ delta*u + iS
 FA =  0.1e1 / x(1) / x(5) ^ 3 * pi ^ (-0.1e1 / 0.2e1) * (-0.4e1 * x(1) ^ 2 * x(5) ^ 3 * sqrt(pi) + 0.3e1 * x(5) ^ 2 * Ib - 0.8e1 * x(1) ^ 2 * x(4) * x(5) ^ 3 * sqrt(pi) - (2 * Id)) / 0.4e1;
 FB = -(0.8e1 * x(1) ^ 2 * sqrt(pi) - 0.5e1 * x(1) ^ 4 * sqrt(0.2e1) * sqrt(pi) * x(5) ^ 2 - 0.4e1 * Is * x(5) ^ 2 + 0.6e1 * Ia * x(5) * x(1) - 0.8e1 * x(1) ^ 2 * x(3) ^ 2 * sqrt(pi) * x(5) ^ 2 - 0.8e1 * x(5) * x(3) * Ic + 0.8e1 * x(1) ^ 2 * delta * sqrt(pi) * x(5) ^ 2) / x(1) ^ 2 / x(5) ^ 2 * pi ^ (-0.1e1 / 0.2e1) / 0.8e1; 
 FC = -(x(3) * Ib + Ixi) / x(1) ^ 2 / x(5) * pi ^ (-0.1e1 / 0.2e1);
 FD = 0.1e1 / x(1) ^ 2 / x(5) ^ 4 * pi ^ (-0.1e1 / 0.2e1) * (0.4e1 * x(1) ^ 2 * sqrt(pi) - 0.16e2 * x(1) ^ 2 * x(4) ^ 2 * x(5) ^ 4 * sqrt(pi) - x(1) ^ 4 * sqrt(0.2e1) * sqrt(pi) * x(5) ^ 2 - 0.4e1 * Is * x(5) ^ 2 + 0.2e1 * Ia * x(5) * x(1)) / 0.4e1;
 FSigma = -(x(5) ^ 2 * Ib - 0.8e1 * x(1) ^ 2 * x(4) * x(5) ^ 3 * sqrt(pi) - (2 * Id)) / x(1) ^ 2 / x(5) ^ 2 * pi ^ (-0.1e1 / 0.2e1) / 0.2e1; 
 FXi = (0.2e1 * x(1) ^ 2 * x(5) * x(3) * sqrt(pi) + Ic) / x(1) ^ 2 / x(5) * pi ^ (-0.1e1 / 0.2e1); 
 
 RHS =  [FA; FB; FC; FD; FSigma; FXi]; 
 VAx(1) = FA;
 VAx(2) = FB;
 VAx(3) = FC;
 VAx(4) = FD;
 VAx(5) = FSigma;
 VAx(6) = FXi;
 %RHS =  [FA; FB; FC; FD; FSigma; FXi]; 
 
end

