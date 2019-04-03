
L=80;
cutoff=L-20;
wcut = (L-cutoff)/2;  % half-width of the cutoff region.

x=linspace(-L,L,1000);

mincut = 0.99;
cutslope = 0.5;

damping_stepM=mincut+(1-mincut)*((1+tanh(cutslope*( x+L-wcut)))/2);
damping_stepP=mincut+(1-mincut)*((1+tanh(cutslope*(-x+L-wcut)))/2);
damping_step=damping_stepM.*damping_stepP;

figure(4);clf
plot(x,damping_step)
axis([-L L 0.9 1+(1-mincut)])

