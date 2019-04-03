function RHS = LLE_steady_state_G4(x)
 
 global delta h sigma_X Ein

 warning('ON','MATLAB:integral:MinStepSize');
 
 FIa = @(t,x) 0.2e1 .* Ein .* exp(-(-t + x(6)) .^ 2 ./ x(1) .^ 2 ./ 0.2e1) .* (sin(x(2) + x(3) .* t - x(3) .* x(6)) .* x(1) .^ 2 + cos(x(2) + x(3) .* t - x(3) .* x(6)) .* t .^ 2 - 0.2e1 .* cos(x(2) + x(3) .* t - x(3) .* x(6)) .* t .* x(6) + cos(x(2) + x(3) .* t - x(3) .* x(6)) .* x(6) .^ 2) ./ x(1) .^ 2;
 Ia = integral(@(t)FIa(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);
 
 FIb = @(t,x) 0.2e1 .* Ein .* x(1) .* exp(-(-t + x(6)) .^ 2 ./ x(1) .^ 2 ./ 0.2e1) .* cos(x(2) + x(3) .* t - x(3) .* x(6));
 Ib = integral(@(t)FIb(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);

 FIc = @(t,x) -0.2e1 .* Ein .* (-t + x(6)) .* x(1) .* exp(-(-t + x(6)) .^ 2 ./ x(1) .^ 2 ./ 0.2e1) .* cos(x(2) + x(3) .* t - x(3) .* x(6));
 Ic = integral(@(t)FIc(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13);

 FIxi = @(t,x) -0.2e1 .* Ein .* exp(-(-t + x(6)) .^ 2 ./ x(1) .^ 2 ./ 0.2e1) .* (-sin(x(2) + x(3) .* t - x(3) .* x(6)) .* t + sin(x(2) + x(3) .* t - x(3) .* x(6)) .* x(6) + x(1) .^ 2 .* cos(x(2) + x(3) .* t - x(3) .* x(6)) .* x(3)) ./ x(1);
 Ixi = integral(@(t)FIxi(t,x),-inf,inf,'AbsTol',1e-15,'RelTol',1e-13); 

 warning('ON','MATLAB:integral:MinStepSize');

 FA = 0.1e1 ./ x(1) .^ 2 .* pi .^ (-0.1e1 ./ 0.2e1) .* (-0.2e1 .* x(1) .^ 3 .* sqrt(pi) + Ib) ./ 0.3e1;
 FB =(0.176e3 .* (x(1) .^ 5) .* (sigma_X .^ 10) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) + 0.128e3 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 2) .* (x(6) .^ 4) + 0.128e3 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 8) .* (x(6) .^ 4) + 0.320e3 .* (x(1) .^ 13) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) .* (sigma_X .^ 2) - 0.256e3 .* (x(1) .^ 15) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* x(6) .* h .^ 2 + 0.320e3 .* (x(1) .^ 13) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) .* (sigma_X .^ 2) + 0.128e3 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 4) .* (sigma_X .^ 2) + 0.48e2 .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) .* (sigma_X .^ 12) + 0.48e2 .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 12) .* (x(6) .^ 2) + 0.384e3 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 4) .* (x(6) .^ 4) + 0.384e3 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 6) .* (x(6) .^ 4) + 0.80e2 .* (x(1) .^ 9) .* (sigma_X .^ 6) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) + 0.176e3 .* (x(1) .^ 7) .* (sigma_X .^ 8) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) + 0.176e3 .* (x(1) .^ 7) .* (sigma_X .^ 8) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) + 0.176e3 .* (x(1) .^ 5) .* (sigma_X .^ 10) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) + 0.224e3 .* (x(1) .^ 11) .* (sigma_X .^ 4) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) + 0.224e3 .* (x(1) .^ 11) .* (sigma_X .^ 4) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) + 0.80e2 .* (x(1) .^ 9) .* (sigma_X .^ 6) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) + 0.384e3 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 4) .* (sigma_X .^ 4) + 0.384e3 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 4) .* (sigma_X .^ 6) + 0.128e3 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 4) .* (sigma_X .^ 8) - 0.512e3 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* x(6) .* (sigma_X .^ 2) - 0.512e3 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* (x(6) .^ 3) .* h .^ 2 .* (sigma_X .^ 2) - 0.1536e4 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* (x(6) .^ 3) .* h .^ 2 .* (sigma_X .^ 4) + 0.256e3 .* (x(1) .^ 7) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 8) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 3) + 0.32e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 10) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 3) + 0.2304e4 .* (x(1) .^ 13) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 4) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* tau0 + 0.2688e4 .* (x(1) .^ 11) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 6) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* tau0 + 0.1536e4 .* (x(1) .^ 9) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 8) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* tau0 + 0.432e3 .* (x(1) .^ 7) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 10) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* tau0 + 0.48e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 12) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* tau0 - 0.512e3 .* (x(1) .^ 13) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 2) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 3) - 0.1024e4 .* (x(1) .^ 11) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 4) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 3) + 0.768e3 .* (x(1) .^ 15) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 2) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* tau0 - 0.768e3 .* (x(1) .^ 15) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 2) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* x(6) - 0.2304e4 .* (x(1) .^ 13) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 4) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* x(6) - 0.2688e4 .* (x(1) .^ 11) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 6) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* x(6) - 0.1536e4 .* (x(1) .^ 9) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 8) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* x(6) - 0.432e3 .* (x(1) .^ 7) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 10) .* x(3) .* sqrt(pi) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* x(6) - 0.48e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 12) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* x(6) + 0.512e3 .* (x(1) .^ 13) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 2) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 3) + 0.1024e4 .* (x(1) .^ 11) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 4) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 3) + 0.768e3 .* (x(1) .^ 9) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 6) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 3) - 0.768e3 .* (x(1) .^ 9) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 6) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 3) - 0.256e3 .* (x(1) .^ 7) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 8) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 3) - 0.32e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 10) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 3) + 0.768e3 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* (x(6) .^ 2) .* h .^ 2 .* (sigma_X .^ 8) + 0.2304e4 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* (x(6) .^ 2) .* h .^ 2 .* (sigma_X .^ 4) + 0.2304e4 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* (x(6) .^ 2) .* h .^ 2 .* (sigma_X .^ 6) - 0.512e3 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* (x(6) .^ 3) .* h .^ 2 .* (sigma_X .^ 8) + 0.768e3 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* (x(6) .^ 2) .* h .^ 2 .* (sigma_X .^ 2) - 0.1536e4 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* (x(6) .^ 3) .* h .^ 2 .* (sigma_X .^ 6) - 0.96e2 .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 .* (sigma_X .^ 12) .* x(6) - 0.640e3 .* (x(1) .^ 13) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 .* (sigma_X .^ 2) .* x(6) - 0.352e3 .* (x(1) .^ 7) .* (sigma_X .^ 8) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* x(6) .* h .^ 2 - 0.352e3 .* (x(1) .^ 5) .* (sigma_X .^ 10) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* x(6) .* h .^ 2 - 0.160e3 .* (x(1) .^ 9) .* (sigma_X .^ 6) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* x(6) .* h .^ 2 - 0.512e3 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* x(6) .* (sigma_X .^ 8) - 0.448e3 .* (x(1) .^ 11) .* (sigma_X .^ 4) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* x(6) .* h .^ 2 - 0.1536e4 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* x(6) .* (sigma_X .^ 4) - 0.1536e4 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* x(6) .* (sigma_X .^ 6) - 0.768e3 .* (x(1) .^ 7) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 8) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 2) .* tau0 - 0.96e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 10) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 2) .* tau0 + 0.1536e4 .* (x(1) .^ 13) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 2) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* x(6) + 0.3072e4 .* (x(1) .^ 11) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 4) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* x(6) + 0.2304e4 .* (x(1) .^ 9) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 6) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* x(6) + 0.768e3 .* (x(1) .^ 7) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 8) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* x(6) + 0.96e2 .* (x(1) .^ 5) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 10) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (tau0 .^ 2) .* x(6) - 0.1536e4 .* (x(1) .^ 13) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 2) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 2) .* tau0 - 0.3072e4 .* (x(1) .^ 11) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 4) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 2) .* tau0 - 0.2304e4 .* (x(1) .^ 9) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (sigma_X .^ 6) .* x(3) .* h .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (x(6) .^ 2) .* tau0 + 0.312e3 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 10) + 0.40e2 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 12) + 0.128e3 .* (x(1) .^ 15) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (x(6) .^ 2) + 0.984e3 .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 8) + 0.192e3 .* (x(1) .^ 17) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* sigma_X + 0.192e3 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sigma_X .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 14) + 0.960e3 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 3) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 12) + 0.2016e4 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 5) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 10) + 0.2304e4 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 7) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 8) + 0.1548e4 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 9) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 6) + 0.612e3 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 11) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 4) + 0.132e3 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 13) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (x(1) .^ 2) + 0.12e2 .* (x(1) .^ 3) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 15) + 0.1548e4 .* (x(1) .^ 9) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 9) + 0.612e3 .* (x(1) .^ 7) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 11) + 0.132e3 .* (x(1) .^ 5) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 13) + 0.2016e4 .* (x(1) .^ 13) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 5) + 0.2304e4 .* (x(1) .^ 11) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 7) + 0.5e1 .* (x(1) .^ 5) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 15) .* sqrt(0.2e1) + 0.192e3 .* (x(1) .^ 17) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* sigma_X + 0.960e3 .* (x(1) .^ 15) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .^ 2 .* (sigma_X .^ 3) + 0.645e3 .* (x(1) .^ 11) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 9) .* sqrt(0.2e1) + 0.255e3 .* (x(1) .^ 9) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 11) .* sqrt(0.2e1) + 0.55e2 .* (x(1) .^ 7) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 13) .* sqrt(0.2e1) + 0.840e3 .* (x(1) .^ 15) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 5) .* sqrt(0.2e1) + 0.960e3 .* (x(1) .^ 13) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 7) .* sqrt(0.2e1) + 0.80e2 .* (x(1) .^ 19) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* sigma_X .* sqrt(0.2e1) + 0.400e3 .* (x(1) .^ 17) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 3) .* sqrt(0.2e1) + 0.1548e4 .* (x(1) .^ 9) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 9) + 0.612e3 .* (x(1) .^ 7) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 11) + 0.132e3 .* (x(1) .^ 5) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 13) + 0.12e2 .* (x(1) .^ 3) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 15) + 0.960e3 .* (x(1) .^ 15) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 3) + 0.2016e4 .* (x(1) .^ 13) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 5) + 0.2304e4 .* (x(1) .^ 11) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* delta .* (sigma_X .^ 7) + 0.1608e4 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 6) + 0.672e3 .* (x(1) .^ 15) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 2) + 0.1440e4 .* (x(1) .^ 13) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 4) + 0.128e3 .* (x(1) .^ 15) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) - 0.320e3 .* (x(1) .^ 13) .* Ia .* (sigma_X .^ 3) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.672e3 .* (x(1) .^ 11) .* Ia .* (sigma_X .^ 5) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.768e3 .* (x(1) .^ 9) .* Ia .* (sigma_X .^ 7) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.516e3 .* (x(1) .^ 7) .* Ia .* (sigma_X .^ 9) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.204e3 .* (x(1) .^ 5) .* Ia .* (sigma_X .^ 11) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.44e2 .* (x(1) .^ 3) .* Ia .* (sigma_X .^ 13) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.4e1 .* x(1) .* Ia .* (sigma_X .^ 15) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.128e3 .* (x(1) .^ 17) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 - 0.258e3 .* (x(1) .^ 7) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 9) - 0.32e2 .* (x(1) .^ 15) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* sigma_X - 0.160e3 .* (x(1) .^ 13) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 3) - 0.336e3 .* (x(1) .^ 11) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 5) - 0.384e3 .* (x(1) .^ 9) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 7) - 0.102e3 .* (x(1) .^ 5) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 11) - 0.22e2 .* (x(1) .^ 3) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 13) - 0.2e1 .* x(1) .* sqrt(pi) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 15) + 0.12e2 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* (sigma_X .^ 15) .* x(3) .* Ic .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.64e2 .* (x(1) .^ 15) .* Ia .* sigma_X .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2))) .* pi .^ (-0.1e1 ./ 0.2e1) .* ((2 .* x(1) .^ 2 + sigma_X .^ 2) .^ (-0.9e1 ./ 0.2e1)) .* ((x(1) .^ 2 + sigma_X .^ 2) .^ (-0.7e1 ./ 0.2e1)) ./ (x(1) .^ 3) ./ sigma_X ./ 0.12e2;
 FC = -(Ixi .* (sigma_X .^ 11) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.32e2 .* (x(1) .^ 11) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.28e2 .* (x(1) .^ 5) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 8) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.4e1 .* (x(1) .^ 3) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 10) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.72e2 .* (x(1) .^ 7) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 6) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.80e2 .* (x(1) .^ 9) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 4) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.24e2 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 .* (sigma_X .^ 6) + 0.16e2 .* (x(1) .^ 7) .* sqrt(pi) .* (x(6) .^ 3) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 2) + 0.32e2 .* (x(1) .^ 5) .* sqrt(pi) .* (x(6) .^ 3) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 4) - 0.8e1 .* (x(1) .^ 7) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 4) - 0.8e1 .* (sigma_X .^ 8) .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* x(6) - 0.16e2 .* (sigma_X .^ 6) .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) + 0.8e1 .* (sigma_X .^ 8) .* (x(1) .^ 3) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 + 0.16e2 .* (sigma_X .^ 6) .* (x(1) .^ 3) .* sqrt(pi) .* (x(6) .^ 3) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 - 0.24e2 .* (sigma_X .^ 2) .* (x(1) .^ 9) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* h .^ 2 + 0.24e2 .* (sigma_X .^ 2) .* (x(1) .^ 9) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 - 0.24e2 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (sigma_X .^ 6) .* x(6) + 0.8e1 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* h .^ 2 .* (sigma_X .^ 4) - 0.16e2 .* (x(1) .^ 7) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* (sigma_X .^ 2) - 0.32e2 .* (x(1) .^ 5) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 3) .* (sigma_X .^ 4) - 0.16e2 .* (x(1) .^ 11) .* sqrt(pi) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* h .^ 2 + 0.16e2 .* (x(1) .^ 11) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 + 0.8e1 .* (x(1) .^ 10) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* x(3) .* sigma_X .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* Ib + 0.28e2 .* (x(1) .^ 8) .* Ib .* x(3) .* (sigma_X .^ 3) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.38e2 .* (x(1) .^ 6) .* Ib .* x(3) .* (sigma_X .^ 5) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.25e2 .* (x(1) .^ 4) .* Ib .* x(3) .* (sigma_X .^ 7) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.8e1 .* (x(1) .^ 2) .* Ib .* x(3) .* (sigma_X .^ 9) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.48e2 .* (x(1) .^ 5) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 6) .* (tau0 .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.48e2 .* (x(1) .^ 5) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (x(6) .^ 2) .* h .* (sigma_X .^ 6) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.48e2 .* (x(1) .^ 7) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) .* (sigma_X .^ 2) + 0.96e2 .* (x(1) .^ 5) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) .* (sigma_X .^ 4) - 0.48e2 .* (x(1) .^ 7) .* sqrt(pi) .* (x(6) .^ 2) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* tau0 .* h .^ 2 .* (sigma_X .^ 2) - 0.96e2 .* (x(1) .^ 5) .* sqrt(pi) .* (x(6) .^ 2) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 .* (sigma_X .^ 4) + 0.48e2 .* (sigma_X .^ 6) .* (x(1) .^ 3) .* sqrt(pi) .* x(6) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* (tau0 .^ 2) - 0.48e2 .* (sigma_X .^ 6) .* (x(1) .^ 3) .* sqrt(pi) .* (x(6) .^ 2) .* exp(-(2 .* (x(6) - tau0) .^ 2 ./ (2 .* x(1) .^ 2 + sigma_X .^ 2))) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) .* h .^ 2 .* tau0 + 0.8e1 .* (x(1) .^ 3) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (x(6) .^ 2) .* h .* (sigma_X .^ 8) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.8e1 .* (x(1) .^ 3) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 8) .* (tau0 .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.64e2 .* (x(1) .^ 9) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (x(6) .^ 2) .* h .* (sigma_X .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.64e2 .* (x(1) .^ 9) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 2) .* (tau0 .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.96e2 .* (x(1) .^ 7) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* h .* (sigma_X .^ 4) .* (tau0 .^ 2) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.96e2 .* (x(1) .^ 7) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* (x(6) .^ 2) .* h .* (sigma_X .^ 4) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) + 0.25e2 .* (x(1) .^ 4) .* Ixi .* (sigma_X .^ 7) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.8e1 .* (x(1) .^ 2) .* Ixi .* (sigma_X .^ 9) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + Ib .* x(3) .* (sigma_X .^ 11) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.8e1 .* (x(1) .^ 10) .* Ixi .* sigma_X .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.28e2 .* (x(1) .^ 8) .* Ixi .* (sigma_X .^ 3) .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) + 0.38e2 .* (x(1) .^ 6) .* (sigma_X .^ 5) .* Ixi .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt((x(1) .^ 2 + sigma_X .^ 2)) - 0.192e3 .* (x(1) .^ 7) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* x(6) .* h .* (sigma_X .^ 4) .* tau0 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.96e2 .* (x(1) .^ 5) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* x(6) .* h .* (sigma_X .^ 6) .* tau0 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.16e2 .* (x(1) .^ 3) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* x(6) .* h .* (sigma_X .^ 8) .* tau0 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi) - 0.128e3 .* (x(1) .^ 9) .* x(3) .* exp(-((x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2))) .* x(6) .* h .* (sigma_X .^ 2) .* tau0 .* sqrt((2 .* x(1) .^ 2 + sigma_X .^ 2)) .* sqrt(pi)) .* ((x(1) .^ 2 + sigma_X .^ 2) .^ (-0.5e1 ./ 0.2e1)) ./ sigma_X .* pi .^ (-0.1e1 ./ 0.2e1) ./ (x(1) .^ 3) .* ((2 .* x(1) .^ 2 + sigma_X .^ 2) .^ (-0.7e1 ./ 0.2e1));
 FXi = (0.2e1 .* sqrt(pi) .* sqrt(x(1) .^ 2 + sigma_X .^ 2) .* x(3) .* x(1) .^ 5 + 0.2e1 .* x(1) .^ 3 .* sqrt(pi) .* sqrt(x(1) .^ 2 + sigma_X .^ 2) .* x(3) .* sigma_X .^ 2 + 0.4e1 .* x(1) .^ 3 .* sigma_X .* sqrt(pi) .* exp(-(x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2)) .* h .* tau0 - 0.4e1 .* x(1) .^ 3 .* sigma_X .* sqrt(pi) .* exp(-(x(6) - tau0) .^ 2 ./ (x(1) .^ 2 + sigma_X .^ 2)) .* h .* x(6) + Ic .* sqrt(x(1) .^ 2 + sigma_X .^ 2) .* x(1) .^ 2 + Ic .* sqrt(x(1) .^ 2 + sigma_X .^ 2) .* sigma_X .^ 2) .* (x(1) .^ 2 + sigma_X .^ 2) .^ (-0.3e1 ./ 0.2e1) .* pi .^ (-0.1e1 ./ 0.2e1) ./ x(1) .^ 3;
 
 RHS =  [FA; FB; FC; FXi]; 
 
end