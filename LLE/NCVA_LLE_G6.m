close all; clear all;

global allt allx 
global delta h sigma_X Ein tau0 uS
config; 
%load('/Users/jules/Documents/Jules/BEC Research /Temporal Tweezing/Code/LLE_PDEv2/SkinnyTweeze/SkinnyLLE.mat')
%load('/Users/jules/Documents/Jules/BEC Research /Temporal Tweezing/Code/LLE_NCVA/RegularTweeze/RegularLLE.mat')
load('/Users/jules/Documents/Jules/BEC Research /Temporal Tweezing/Code/LLE_NCVA/FatTweeze/FatLLE.mat')

%load('RegularLLE.mat')
%load('FatLLE.mat')
delta
sigma_X
h
%a = 50;
%b = 100;
%r = (b-a).*rand(1000,1) + a;
% [a, b, c, d, s, xi]
% a_random = (10-0).*rand(100,1) + 0;
% b_random = (2+2).*rand(100,1) - 2;
% d_random = (2+2).*rand(100,1) -2;
% s_random = (5-0).*rand(100,1) + 0 ;
%for i = 2:length(a_random);
i = 1;
    x0 =    [ 3.3366; 0.3108; 0.0000; 0.8486; 0.6836; 0.0000];
    x0 = fsolve(@LLE_steady_state_G6B, x0);
    [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_G6B(solp_),1e-15*[1,1],sol_parms);
    %[uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,eta,dx,xi0),1e-15*[1,1],sol_parms);
    %[solp, resnorm,Ein residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[],options);
    %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
    uVA(:,i) = (solp(1)).*exp(-0.5*(((xi0-solp(6))/solp(5)).^2)).*exp(1i*(solp(4).*(xi0-solp(6)).^2 + solp(3).*(xi0-solp(6)) + solp(2)));
    %thetauVA = atan2(imag(uVA(:,i)), real(uVA(:,i)));
    %uVA(:,i) = uVA(:,i).*exp(-1i*thetauVA);
    params(:,i) = [solp(1); solp(2); solp(3); solp(4); solp(5); solp(6)];
    %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
    bVA(i) = uVA(151,i).*conj(uVA(151,i));
    fRHS(:,i) = LLE_steady_state_G6B(solp');
    test = sum(fRHS);
    error(i) = ierr;
%end 
p0 = x0;

save('FatNCVA','params','p0','uVA','delta','h','sigma_X','xi0','tau0','uS')