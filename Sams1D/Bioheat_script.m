%This is a run script for the Bioheat1D.m file. First load and parse the
%power
cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/010/laser_log
load 'power_log.txt';
[P,unique_P,delta_P]=power_parser(power_log);
%P=[1 -3; 5 0;15 3; 20 6; 25 9; 30 12; 35 15];

%Define the domain and scaling
mod_point.x=51;
mod_point.y=51;
mod_point.z=1;

matrix.x=256;
matrix.y=256;
matrix.z=1;

scaling.x=2;
scaling.y=2;
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

%Define the source; source.n must be greater than 1 to make sense. Odd is
%better than even
source.n=5;
source.length=0.035;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

%Run the Bioheat model with the unique powers
tic
[tmap_unique]=Bioheat1D_zeroC(P,dom,source);
toc

%tmap_unique=tmap_unique-273.15-(1277.6-37); %Stupid correction for unknown bias.
tmap_unique=tmap_unique+37;
tmap_unique(:,:,:,1)=37;
%[tmap]=Build_tmap_history(tmap_unique,delta_P);

clear FOV scaling matrix mod_point power_log delta_P mod_point unique_P source mod_pix MRTI_pix dom



figure(1);imagesc(tmap_unique(:,:,1,2));


%Make the full temperature history with the unique tmaps
%[tmap]=Build_tmap_history(tmap_unique,P,unique_P);