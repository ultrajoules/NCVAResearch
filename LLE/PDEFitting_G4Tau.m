%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit PDE to ansatz 2-p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    %Compare parameters a, b, c, xi from ODE solution to PDE -
    %this is done by projecting the PDE solution to ODE ansatz by
    %fitting at select time intervals...
    %1 = a;  2 = b;  3= c; 4 = d; 5 = sigma; 6 = xi;
    pp = zeros(length(allt), 4);
    gaussian4_func = inline('[(1/sqrt(2*pi*ppde(1)^2)).*exp(-0.5*((x-ppde(4))/(ppde(1))).^2).*cos((ppde(2) + ppde(3).*(x-ppde(4))));(1/sqrt(2*pi*ppde(1)^2)).*exp(-0.5*((x-ppde(4))/(ppde(1))).^2).*cos((ppde(2) + ppde(3).*(x-ppde(4))))]','ppde','x');
    xLS = xi0;
    p0 = [sqrt(5), 1, 1,  tau0 ];
    pp = p0;
    options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'MaxFunEvals', 100000, 'MaxIter', 100000);
    %options.MaxIter = 1000;
    %options = optimset('Display','off','TolFun',1e-10);
    [ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(gaussian4_func,p0,xLS,[real(u);imag(u)],[],[],options);
    uppde = (1/sqrt(2*pi*ppde(1)^2)).*exp(-0.5*((xi0-ppde(4))/(ppde(1))).^2).*exp(1i*(ppde(2) + ppde(3).*(xi0-ppde(4))));
%     for jj = 1:saveframes+1
%         yLS(:) = allU(jj,:);
%         [ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(sech_func,p0,xLS,[real(yLS);imag(yLS)],[],[],options);
%         pp(jj,:) = ppde;
%         p0 = ppde;
%     end
    toc;