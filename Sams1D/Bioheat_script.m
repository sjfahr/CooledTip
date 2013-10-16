%This is a run script for the Bioheat1D.m file. First load and parse the
%power
cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/010/laser_log
load 'power_log.txt';
[P,unique_P,delta_P]=power_parser(power_log);

%Define the domain and scaling
mod_point.x=51;
mod_point.y=51;
mod_point.z=1;

matrix.x=256;
matrix.y=256;
matrix.z=1;

scaling.x=1;
scaling.y=1;
scaling.z=1;

FOV.x=0.25; %change for image/dataset
FOV.y=0.25;
FOV.z=0.007; 

%Build the domain
[dom,MRTI_pix,mod_pix]=modeled_domain(FOV,matrix,scaling,mod_point);
%Display the domain
dom
MRTI_pix
mod_pix

%Define the source
source.n=10;
source.length=0.01;

%Run the Bioheat model with the unique powers
[tmap_unique]=Bioheat1D(P,dom,source);

tmap_unique=tmap_unique-273.15;

[tmap]=Build_tmap_history(tmap_unique,delta_P);

clear FOV scaling matrix mod_point power_log delta_P mod_point unique_P tmap_unique source mod_pix MRTI_pix dom


%Make the full temperature history with the unique tmaps
%[tmap]=Build_tmap_history(tmap_unique,P,unique_P);