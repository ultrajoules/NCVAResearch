%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NLSE variable definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EinVec = 2; 
hVec = 2; 
sigma_XVec = 0.5;
tau0 = 0; 

damping_step=1;
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
Tend =5; 
dt = 0.0001;
maxTiterations = Tend/dt;
xi0 = x';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
mu = 0;
alpha = -1; 
maxit=10000;                          %nsoli params
maxitl=10000;
etamax=0.9;
lmeth=3;
restart_limit=10000;
sol_parms=[maxit,maxitl,etamax,lmeth,restart_limit];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

