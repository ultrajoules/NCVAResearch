%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit PDE to ansatz 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    %Compare parameters a, b, c, xi from ODE solution to PDE -
    %this is done by projecting the PDE solution to ODE ansatz by
    %fitting at select time intervals...
    %1 = a;  2 = b;  3= c; 4 = xi;
    pp = zeros(length(allt), 4);
    func2 = inline('real(ppde(1)).*sech(real(ppde(1)).*(x-real(ppde(4)))).*exp(1i*(real(ppde(3))*(x - real(ppde(4))) + real(ppde(2))))','ppde','x');
    sech_func = inline('[ppde(1).*sech(ppde(1).*(x-ppde(4))).*cos(ppde(3)*(x-ppde(4))+ppde(2));ppde(1).*sech(ppde(1).*(x-ppde(4))).*sin(ppde(3)*(x-ppde(4))+ppde(2))]','ppde','x');
    sechexp_func = inline('[ppde(1).*sech(ppde(1).*(x-ppde(4))).*cos(ppde(3)*(x-ppde(4))+ppde(2).*(x-ppde(4)).^2);ppde(1).*sech(ppde(1).*(x-ppde(4))).*sin(ppde(3)*(x-ppde(4))+ppde(2).*(x-ppde(4)).^2)]','ppde','x');
    xLS = xi0;
    p0 = [max(real(u)), 1, 0, tau0];
    pp = p0;
    options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'MaxFunEvals', 100000, 'MaxIter', 100000);
    %options.MaxIter = 1000;
    %options = optimset('Display','off','TolFun',1e-10);
    [ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(sechexp_func,p0,xLS,[real(u);imag(u)],[],[],options);
%     for jj = 1:saveframes+1
%         yLS(:) = allU(jj,:);
%         [ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(sech_func,p0,xLS,[real(yLS);imag(yLS)],[],[],options);
%         pp(jj,:) = ppde;
%         p0 = ppde;
%     end
    toc;