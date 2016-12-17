
clc; clear;
N_sat=21;   %number of sats in a plane
N_planes=9; %number of orbital planes
h=542;      %high of the sats in km
I=72;       %inclination of the orbits in degrees
phase=210/(N_planes-1); %Angle between planes in the equator
At=0.1;       %time step in minutes
T=48;       %time to simulate in hours
t=0:At:T*60;    %array of time in minutes
e_min=7.5;       %minumum elevation in degrees
lambda=57.5; %latitudes to simulate in degrees
lat=length(lambda);
mu=0;           %longitudes to simulate in degrees
long=length(mu);
l=length(t);

Xs=Orbital_position(N_sat,N_planes,h,I,phase,t);
L=zeros(lat,l);

for i=1:lat
    Xg=Ground_position(lambda(i),mu,t);
    [L(i,:)]=Links(Xs,Xg,e_min);
end

%%
figure; plot(t/60,L); legend 0 30 60 90; xlabel('time (h)'); ylabel('links');

