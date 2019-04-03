%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NLSE variable definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ein = [1.61:-0.01:0]; %sqrt(.960)*sqrt(1.2*0.1*0.1/0.146^3);
EinVec = 2; %sqrt(.960)*sqrt(1.2*0.1*0.1/0.146^3); %2; %[1.61:-0.01:0]; [0.5:0.5:10]; %[0:0.5:10];
%hVec = 2; %[0:0.5:10];%[0:0.5:5];
sigma_XVec = 3; %0.5; %[0.5:0.5:10]; %[0:0.5:10];
%sigma_X = 9.0/(2*sqrt(2*log(2))); %2.3/sqrt(2); %2;  %.90/(2*ln(2)2;
tauSpeed = @(t) -0.01*t;
tau0 = 0; %tauSpeed(0);
delta = 0.41/0.146; 
%deltaVec = [1:0.5:5];
eta = -1;
onehalf = 1.0;  %FOr LLE should be 1
g = -1; %Bright/focusing/attractive = -1,   Dark/defocusing/repulsive = +1
E = @(x, tau0, sigma_X, h) h.*exp(-0.5*((x-tau0)/sigma_X).^2);
dE = @(x, tau0, sigma_X, h) -(h.*(x-tau0)./sigma_X^2).*exp(-0.5*((x-tau0)/sigma_X).^2);
ddE = @(x, tau0, sigma_X, h) -h.*exp(-0.5*((x-tau0)/sigma_X).^2)./sigma_X^2 + (h.*(x-tau0).^2./sigma_X^4).*exp(-0.5*((x-tau0)/sigma_X).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define discetization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax = 40;
dx = 0.05; 
oodx2=1/dx^2;
oodx = 0.5/dx;
oos = 1/6;
dxx = 1/dx^2;
ohdxx = onehalf*dxx;
onehalfoodx2=onehalf/(dx*dx);
x = [-xmax:dx:xmax];
N = length(x);
t0 = 0;
Tend =10; 
dt = 0.0001;
%dt = 0.0001;
maxTiterations = Tend/dt;
xi0 = x';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hVec = (sqrt(2)*sigma_XVec^2)/max(xi0.*exp(-xi0.^2/(2*sigma_XVec^2))); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define for saving frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveframes = 20;
maxstep = Tend/dt;
stopsave = fix(maxstep/saveframes);

allx = xi0;
allt = t0:dt*stopsave:Tend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define NSOLI params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0; %s0*max(real(uini));
alpha = -1; 
maxit=1000;                          %nsoli params
maxitl=1000;
etamax=0.9;
lmeth=3;
restart_limit=1000;
sol_parms=[maxit,maxitl,etamax,lmeth,restart_limit];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables for ansatz initial conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-p Sech Ansatz %
    uinit = @(x,t,a,w,b,c,d,xi) a.*sech(w*(x-xi)).*exp(1i*(d.*(x-xi).^2 + c.*(x-xi) + b));
  %  fitting_func =  inline(...
  %      '[(p(1)).*sech(p(1)*(x-p(4))).*exp(1i*(p(3).*(x-p(4)) + p(2).*(x-p(4))^2))]'...
  %      ,'p','x');
%4-p SechExp Ansatz
    fitting_func =  inline(...
        '[(p(1)).*sech(p(1)*(x-p(4))).*exp(1i*(p(3).*(x-p(4)) + p(2).*(x-p(4)).^2))]'...
        ,'p','x');
    
%2-p Gaussian function 
    fitting_func2 = inline(...
          '[(1/sqrt(2*pi)).*exp(-0.5*((xi0-solp(2)).^2)).*exp(1i*(solp(1).*(xi0-solp(2))))'...
          ,'p','x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
