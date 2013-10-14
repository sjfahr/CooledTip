%This is a run script for the Bioheat1D.m file.

P=[1 5]; %Power

%Define the domain
mod_point=[51 51 1];

matrix=[256 256 1];

scaling=[1 1 1];

FOV=[0.25 0.25 0.007];  %change for image/dataset

[dom,dom_point,MRTI_pix,mod_pix]=modeled_domain_array(FOV,matrix,scaling,mod_point);

dom
dom_point
MRTI_pix
mod_pix