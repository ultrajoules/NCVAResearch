    i=1;
    x0 = [1;0;0;tau0]; %ppde'; %+randn(4,1); %[a_random(i); b_random(i); c_random(i); xi_random(i)];
    [solp,it_hist,ierr,u_hist]=nsoli(x0,@(solp_)LLE_steady_state_G4Tau(solp_),1e-15*[1,1],sol_parms);
    %[uNSOLI,it_hist,ierr,u_hist]=nsoli(u0,@(uNSOLI_)GPE_nsoli_Stationary_1D(uNSOLI_,N,onehalfoodx2,V,eta,dx,xi0),1e-15*[1,1],sol_parms);
    %[solp, resnorm,Ein residual, exitflag] = lsqnonlin(@LLE_steady_state_XC6,x0,[],[],options);
    %[solp, fval, exitflag, output] = fsolve(@PGGaussian,x0)
    uVA(:,i) = (1/sqrt(2*pi*solp(1)^2)).*exp(-0.5*((xi0-solp(4))/(solp(1))).^2).*exp(1i*(solp(2) + solp(3).*(xi0-solp(4))));
    %thetauVA = atan2(imag(uVA(:,i)), real(uVA(:,i)));
    %uVA(:,i) = uVA(:,i).*exp(-1i*thetauVA);
    params(:,i) = [solp(1); solp(2); solp(3); solp(4)];
    %leastSquares(i) = sum((uVA(:,i).*conj(uVA(:,i)) - uPDE(:).*conj(uPDE(:))).^2);
    bVA(i) = uVA(151,i).*conj(uVA(151,i));
    fRHS(:,i) = LLE_steady_state_G4Tau(solp');
    error(i) = ierr;