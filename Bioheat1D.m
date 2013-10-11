%As a function of power, P, this makes a pre
function [tmap]=Bioheat1D(P,x-dom,y-dom,z-dom,);
R1=0.0015/2; % (m) R1 is the distance from the isotropic laser source point and the edge of the fiber
R2=1; % (m) R2 is the maximum edge of the domain;

mua=500; % 1/m
mus=14000; %1/m
g=0.88; % Unity
k=0.527; % W/(m * K)
w=6; % kg / (m^3 * s)
u0=37+273.15; % K
ua=37+273.15; % K