%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% ORBIT DECAY 

% This routine computes the orbital decay of an spacecraft with time
% using the cowell's method. 
% The differential equations are integrated using "ode45"(Runge-Kutta 4th,
% 5th order)
 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 
clear all;
clc;
close all;

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global mu mmu smu re rs rm omega J2;

%_______Physical Variables____________________________
mu=3.986004418e14; % earth gravitational constant [m^3/s^2]
mmu=4902.800076e9; %moon gravitational constant [m^3/s^2]
smu=132712440040.944e9; %sun gravitational constant [m^3/s^2]

re=6378136.3; %earth radius [m]
rs=696e6; %sun radius [m]
rm=1730e3; %moon radius [m]

omega=7.292115e-5; %earth angular velocity [rad/s]
J2=0.001081874; %earth oblateness gravity coeff []

au=149597870691; %astronomical unit [m]
c=2.99792458e8; %[m/s] speed of light

%____Numerical variables_______
hlim=180; %minimum height for considering the decay [km]
dt=15*60; %timestep of the simulation [s]
t0=0; %initial time of the sim
tf=dt;

%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__Orbital elements______________
a=re+542e3; %[m] semimajor axis
e=0.0001; % eccentricity
i=degtorad(72); %[rad] inclination
w=0; %[rad] argumetn f periapsis
Omega=0; %[rad] RAAN
nu=0; %true anomaly

date=datetime(2017,2,4); %datetime array [Y, M, D]

%__Oblateness parameters_________
znls=0; %zonals [0-18]
tssrls=0; %tesserals [0-18]

%__Drag perturbation inputs______
Cd=2.2; %Drag coefficient (Typical for sats, various ref.)
Adrag=.09; %[m^2] Wet area for the drag
m=4.1; %[kg] Cubesat mass

Bc=Adrag*Cd/m; %Ballistic coeff.

%__SRP perturbation inputs______
Cr=2; %reflectivity ct.
Asrp=.5; %[m^2] Area for the SRP

% Ws (@Dap) = 1361/(1+0.0334*cos(2*pi*Dap/365);
Ws=1361; %W/m^2 approximation for the solar irradiance

srp=Ws/c*Cr*Asrp/m*au^2; %Solar radiation pressure ct.

%% Main Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Computing...');
f=@(q) q*((3+3*q+q^2)/(1+(1+q)^3/2)); %richard Battin's function

aobl=zeros(3,1); 
adrag=aobl; asun=aobl; amoon=aobl; asrp=aobl;


y=cartesian2(a,e,i,Omega,w,nu); %obtain the state vector (r,V)
jdate0=juliandate(date); %compute the julian date

tic
iter=0; tt=[]; xx=[];

h=(norm(y(1:3))-re)*1e-3;                                        
while h>hlim
    jdate=jdate0+t0/86400;
%     GMST=JD2GMST(jdate);
%     
%     %Acceleration due to oblateness
%     aobl=Oblat((t0+tf)/2,y,znls,tssrls,GMST);
%     
%     %Acceleration due to the 3rd body pert.
%     rsun=planetEphemeris(jdate,'Earth','Sun');
      %geocentric pos of the Sun
%     rmoon=planetEphemeris(jdate,'Earth','Moon'); 
      %geocentric pos of the Moon
%     rsun=rsun'; rmoon=rmoon';
%     
%     rs2s=y(1:3)-rsun; %vector sun-sat
%     rm2s=y(1:3)-rmoon; %vector moon-sat
%     
%     qs=dot(y(1:3),(y(1:3)-2*rsun))/dot(rsun,rsun);
%     qm=dot(y(1:3),(y(1:3)-2*rmoon))/dot(rmoon,rmoon);
%     
%     asun=-smu/norm(rs2s)^3*(y(1:3)+f(qs)*rsun);
%     amoon=-mmu/norm(rm2s)^3*(y(1:3)+f(qm)*rmoon);
%     
%     %Acceleration due to the solar radiation presure.
%     usun=rsun/norm(rsun); %unit vector of the pos of the sun
%     us=y(1:3)/norm(y(1:3)); %unit vector of the pos of the sat
%     
%     shadow=random('Normal',0.5,.1); %shadow factor 
%     
%     asrp=srp*rs2s/norm(rs2s)^3;
    
    %Acceleration due to the drag force
    % Compute exospheric temperature [K]
    T = 900 + 2.5*(100-70);
    % Compute effective atmospheric molecular mass [km/K], valid 180<H<500
    M = 27 - 0.012 * (h - 200);
    % Compute atmospheric scale height [km]
    SH = T / M;
    % Compute atmospheric density [kg/m3]
    rho = 6E-10 * exp(-(h - 175) / SH);
    
    atmosV=y(4:6)+omega*[y(2); -y(1); 0]; 
    %relative vel between sat and atmos
    
    adrag=-.5*rho*Bc*atmosV*norm(atmosV);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    acc=aobl+asun+amoon+asrp+adrag;
    
    ode=@(t,y)[y(4);
            y(5); 
            y(6);
            acc(1)-mu/norm(y(1:3))^3*y(1);
            acc(2)-mu/norm(y(1:3))^3*y(2);
            acc(3)-mu/norm(y(1:3))^3*y(3)
            ]; %function to integrate (6eqn)
    
    [t,x]=ode45(ode,[t0 tf], y);
    
    y=x(size(x,1),:);
    %save one value each 12 hours
    if mod(t(length(t)),12*3600)==0  
        tt=[tt; t(length(t))]; xx=[xx; y];
    end
    y=y';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %End of the timestep
    t0=tf;
    tf=tf+dt;
    
    h=(norm(y(1:3))-re)*1e-3;
    
    iter=iter+1;
    if mod(iter,300)==0
        fprintf('Days elapsed: %0.1f \nH=%0.3f km\n',(jdate-jdate0),h);
    end
end
toc
%% Post Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('The satellite decays in %0.0f days',tf/86400);

h=zeros(size(xx,1),1);
for i=1:size(xx,1)
    h(i)=(norm(xx(i,1:3))-re)*1e-3;
end

plot(tt/86400,h);

