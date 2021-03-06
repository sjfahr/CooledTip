%This script is meant to be run in SI units.
%As a function of power, P, this makes a prediction of the temperature
%within the domain, dom. 'dom' is a structure representing the simulated domain with
%three spacial dimensions in meters and the number of points in each
%dimension (a total of 6 fields). The number of points in 'dom' should be
%odd in order for the centroid to be on the origin and a point to be on the
%origin.

%'source' is a structure with the details of the source. The 'n' field is
%the number of sources. The 'length' is the diffusing tip length, commonly
%0.01 m.

function [tmap]=Bioheat1D(P,dom,source,w,k,g,mua,mus,probe_u,robin_co);

%List of space and time details
R1=0.0015/2; % (m) R1 is the distance from the isotropic laser source point and the edge of the fiber
R2=1; % (m) R2 is the maximum edge of the domain;
Npowers=size(P,1)+1; %Returns how many timesteps will be calculated
P(2:Npowers,:)=P(1:(Npowers-1),:); %Inserts the P=0 timestep
P(1,:)=0;
P(:,2)=P(:,2)/source.n;  %Scales the power to the number of source points

%List of constants
% mua=500; % 1/m
% mus=14000; %1/m
% g=0.88; % Unity
% k=0.527; % W/(m * K)
% w=6; % kg / (m^3 * s)
u0=0; % K
ua=0; % K

%Points structure
points.x=linspace(-dom.x/2,dom.x/2,dom.pointx);
points.y=linspace(-dom.y/2,dom.y/2,dom.pointy);
points.z=linspace(-dom.z/2,dom.z/2,dom.pointz);

%Initialize tmap, t_sample, and r

t_sample=zeros(dom.pointx,dom.pointy,dom.pointz,source.n,Npowers);
r=zeros(source.n,1);

%Spatial locations of the sources; My convention is that the long axis of
%the laser is parallel to the y-axis.
%laser=linspace((-source.length/2),(source.length/2),source.n);

%Giant for loop vector for each source, calculate the t_sample; i, ii, iii
%are spatial; j and jj are non-spatial
for i=1:dom.pointx   %Spatial loop for i, ii, iii
    
    for ii=1:dom.pointy
        
        for iii=1:dom.pointz
            
            for j=1:source.n    %Loop for the separate isotropic sources
                r(j)=sqrt(points.x(i)^2+(source.laser(j)-points.y(ii))^2+points.z(iii)^2);  %Distance for each isotropic source
                
                for jj=2:Npowers  %Loop for the unique power settings;.
                    [t_sample(i,ii,iii,j,jj)]=sammodel1D(u0,ua,k,w,P(jj,2),r(j),mua,mus,R1,R2,g);
                    
                end
            end
        end
    end
end

tmap=sum(t_sample,4);  %sum(t_sample(j,jj,1));

end