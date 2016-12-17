function [pert rho] = Perturbation( a,e,i,Bc)
%% Perturbations
% This matlab is used to compte the perturbaons in the classical orbital
% elements due to:
% � J2 (non spherical earth)
% � Atmospheric drag
% � Third body (Moon and Sun)
%Inputs:
%   a: semi-major axis [m]
%   e: eccentricity [-] 
%   i: inclination [rad]
%   m: mass of the satellite [kg]
%   Cd: drag coefficient [-]
%   A: Area of the satellite  [m^2]
%Outputs:
%   Pert: matrix (4x5) with the perturbations causes at the rows (J2,Drag,
%   Moon,Sun) and the orbital elements in the columns (a,e,i,omega,Omega)
%   incP/revolution
%   incv/revolution
%TO DO:
%Les perturbacions degudes a J2 i 3rd body no se si son en pert/dia o no,
%cal mirar-ho per tal de tenir totes les pertrubacions escalades igual
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Previous calculations and PreAllocation
%Physical data
RE=6.378e6; %Earth Radius [m]
u=3.986e14; %GM Earth
J2=1082.6e-5; %J2 coefficient;

P=2*pi*(a^3/u)^.5; % Orbit Period [s]
v=sqrt(u/a); %satellite velocity [m/s]
n=(2*pi/P)*86400;  %number of revolutions/day


% ATMOSPHERIC MODEL
%get the density from MSISE model (extracted as a list from 
%http://omniweb.gsfc.nasa.gov/vitmo/msis_vitmo.html
%    h[km]    rho[kg/m^3]
% filename1='msis_26371.lst'; filename2='msis_2226.lst';
% formatSpec='%f %f';
% size=[2 Inf];
% fileID=fopen(filename1,'r');
% rho_h=fscanf(fileID,formatSpec,size);
% fclose(fileID);
% fileID=fopen(filename2,'r');
% rho_h=[rho_h fscanf(fileID,formatSpec,size)];
% fclose(fileID);
% rho=interp1(rho_h(1,:),rho_h(2,:)',(a-RE)/1000);
    
    H = (a - RE)*1e-3;
    % Compute exospheric temperature [K]
    T = 900 + 2.5*(120-70);
    % Compute effective atmospheric molecular mass [km/K], valid 180<H<500
    M = 27 - 0.012 * (H - 200);
    % Compute atmospheric scale height [km]
    SH = T / M;
    % Compute atmospheric density [kg/m3]
    rho = 6E-10 * exp(-(H - 175) / SH);
    

pert=zeros(4,5);

%% Computation

%J2 pertubation________________ PER REVOLUCIO
j=1;
pert(j,5)=-1.5*n*J2*(RE/a)^2*cos(i)*(1-e^2)^-2;
pert(j,4)=0.75*n*J2*(RE/a)^2*(4-5*(sin(i))^2)*(1-e^2)^-2;

%Atmospheric Drag______________ PER REVOLUCIO
%Perturbations/revolution
j=2;
pert(j,1)=-2*pi*rho*a^2/Bc;
 
% incP=-6*pi^2*a^2*rho/(v*Bc);
% incv=pi*a*rho*v/Bc;

%Moon Perturbation_____________ /DIA --> /REV
j=3;
pert(j,5)=-0.00338*cos(i)/n;
pert(j,4)=0.00169*(4-5*(sin(i))^2)/n;

% Fins aqui son �/dia, 
% Canviem unitats + apliquem periode actual
pert(j,5)=pert(j,5)*P/(24*3600);
pert(j,4)=pert(j,4)*P/(24*3600);

%Sun Perturbation_____________ /DIA --> /REV
j=4;
pert(j,5)=-0.00154*cos(i)/n;
pert(j,4)=0.00077*(4-5*(sin(i))^2)/n;

% Fins aqui son �/dia, 
% Canviem unitats + apliquem periode actual
pert(j,5)=pert(j,5)*P/(24*3600);
pert(j,4)=pert(j,4)*P/(24*3600);


end