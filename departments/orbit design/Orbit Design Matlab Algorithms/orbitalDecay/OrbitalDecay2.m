%  This matlab routine computes the orbital decay of a spacecraft without T,
%  only by taking into account the contributins of the atmospheric drag. 
%
%  Given the satellite initial altitude, the trajectory of the satellite 
%  will be  computed.

% ___References:_____
% [1] Pamrar, R. Satellite Orbital Decay Calclulations. The Australian Space Weather
% [2] An Evaluation of CubeSat Orbital Decay.pdf
% [3]
%%
clc; 
clear all;
close all;
%% 1. Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%___Simulation Parameters____
h=1;   %Time step in terms of period (dt=h*P)

%___Physical Parameters____
u=3.986e14; %GM Earth
RE=6.378e6; %Earth Radius [m]
p=500e3; %Parking orbit height [m]

%Spacecraft phyical data
m0=4; %Inicial spacecraft mass [kg]
dim=[10 10 30]*1e-2;  %dimensions of the cubesats (in cm) [width height long]
cd=2.2;   %drag coefficient---> [1]

%% 2. Previous Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Previous calculations...');

Rp=RE+p;
S=1/2*(dim(1)*dim(2)+dim(1)*dim(3)+dim(2)*dim(3)); %averaged surface model-> extr. [1]
S=0.1;
%get the density from MSISE model (extracted as a list from 
%http://omniweb.gsfc.nasa.gov/vitmo/msis_vitmo.html
%    h[km]    rho[kg/m^3]
filename1='msis_26371.lst'; filename2='msis_2226.lst';
formatSpec='%f %f';
size=[2 Inf];
fileID=fopen(filename1,'r');
rho_h=fscanf(fileID,formatSpec,size);
fclose(fileID);
fileID=fopen(filename2,'r');
rho_h=[rho_h fscanf(fileID,formatSpec,size)];
fclose(fileID);


rho_h(2,:)=rho_h(2,:)*1e3; %pass from g/cm^3 to kg/m^3

%___Preallocation_________
r=[]; P=[]; t=[];

%% 3. initial Map
t(1)=0; 
r(1)=Rp; %Satellite Orbital radius [m]
P(1)=2*pi*(r(1)^3/u)^.5; % Orbit Period [s]
dt=h*P(1);
 
%% 4. Solve the time steps
display('Solving..');

%Iterate satellite with time until it achieves the ground
i=1;
 while r(i)>RE
     i=i+1;
     
     p=r(i-1)-RE; %Sat height
     rho=interp1qr(rho_h(1,:),rho_h(2,:)',p/1000);  % compute the density
     
     if isnan(rho) break; end
     
     %dP=-6*pi^2*rho*r(i-1)^(5/2)/sqrt(u)*(cd*S/m0)*dt;  % Reduction in the period/rev
     dP=-3*pi*rho*r(i-1)*(cd*S/m0)*dt;  % Reduction in the period/rev
     
     %Compute the variation on the orbital period:
     P(i)=P(i-1)+dP;
     r(i)=(P(i)^2/(4*pi^2)*u)^(1/3); % Orbit Period [s]

     t(i)=t(i-1)+dt;
     dt=h*P(i);
     
     %show one timestep each 4000 iterations
    if mod(i,4000)==0
        fprintf('\nDías t=%g  h=%g km\n', t(i)/(86400),(r(i)-RE)/1000);
    end
    
    if r(i)>r(i-1)||r(i)<RE
         r(i)=NaN;
         break;
     end
 end
 
%% 5. PostProcess
fprintf('\nGround achieved after %g days\n', t(i)/86400);
plot(t/86400,(r-RE)/1e3);
grid;
xlabel('time [days]');ylabel('altitude [km]');
tit=['Satellite Orbital decay (m=' num2str(m0) ' kg, S_a_v_g=' num2str(S) ' m^2)']; 
title(tit);