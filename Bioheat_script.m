%This is a run script for the Bioheat1D.m file.

P=[1 5]; %Power

%Define the domain
dom.x=0.2;
dom.y=dom.x;
dom.z=0;
dom.pointx=51;
dom.pointy=dom.pointx;
dom.pointz=1;

%Define the source
source.n=5;
source.length=0.01;

[tmap]=Bioheat1D(P,dom,source);