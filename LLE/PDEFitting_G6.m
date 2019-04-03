%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit PDE to ansatz 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    %Compare parameters a, b, c, xi from ODE solution to PDE -
    %this is done by projecting the PDE solution to ODE ansatz by
    %fitting at select time intervals...
    %1 = a;  2 = b;  3= c; 4 = d; 5 = sigma; 6 = xi;
    pp = zeros(length(allt), 6);
    gaussian6_abs = inline('(pde(1).^2).*exp(-0.5*((x-pde(3))/pde(2)).^2).^2','pde','x');
    gaussian6_func = inline('[ppde(1).*exp(-0.5*((x-ppde(6))/ppde(5)).^2).*cos(ppde(2) + ppde(3)*(x-ppde(6))+ppde(2).*(x-ppde(6)).^2);ppde(1).*exp(-0.5*((x-ppde(6))/ppde(5)).^2).*sin(ppde(2) + ppde(3)*(x-ppde(6))+ppde(2).*(x-ppde(6)).^2)]','ppde','x');
    xLS = xi0;
    %p0 = [sqrt(8), 5, 0, 0, 0.75,  tau0];
    %gaussian6_func2(ppde, x, a, s, xi)  
    p0 = [max(real(u)), 0.75,  tau0];
    options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'MaxFunEvals', 100000, 'MaxIter', 100000);
    pp = p0;
    [pde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(gaussian6_abs,p0,xLS,[u.*conj(u)],[],[],options);
    a = pde(1);
    s = pde(2);
    xi = pde(3);
    p0 = [0, 0, 0];
    g6_func2 = @(ppde, x)gaussian6_func2(ppde,x,a,s,xi);
    pp = p0;
    [ppde2,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(g6_func2,p0,xLS,[real(u);imag(u)],[],[],options);
    ppde = [a; ppde2(1); ppde2(2); ppde2(3); s; xi];
    %options.MaxIter = 1000;
    %options = optimset('Display','off','TolFun',1e-10);
    %[ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(gaussian6_func,p0,xLS,[real(u);imag(u)],[],[],options);
    uPDE = (ppde(1)).*exp(-0.5*(((xi0-ppde(6))/ppde(5)).^2)).*exp(1i*(ppde(4).*(xi0-ppde(6)).^2 + ppde(3).*(xi0-ppde(6)) + ppde(2)));

%     for jj = 1:saveframes+1
%         yLS(:) = allU(jj,:);
%         [ppde,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(sech_func,p0,xLS,[real(yLS);imag(yLS)],[],[],options);
%         pp(jj,:) = ppde;
%         p0 = ppde;
%     end
    toc;
    

