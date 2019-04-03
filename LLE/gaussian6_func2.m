    function y = gaussian6_func2(ppde, x, a, s, xi)  
     y = [a.*exp(-0.5*((x-xi)/s).^2).*cos(ppde(1) + ppde(2)*(x-xi)+ppde(3).*(x-xi).^2);a.*exp(-0.5*((x-xi)/s).^2).*sin(ppde(1) + ppde(2)*(x-xi)+ppde(3).*(x-xi).^2)]; 
    end