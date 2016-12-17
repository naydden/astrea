%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% ORBIT PLOTTER - Influential phenomena computation

% This matlab routine plots the position of the satellites of a Walker
% Delta / SemiWalker Delta or other generated configurations.

clc, clear, close all;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Data Input
%add the subfolder programs to the current path
addpath(genpath('./auxiliar'));

%_____Physical input__________

global h, h = 542;  % [km]
degreegen = 225;    % Walker Delta 360�, SemiWalker Delta 180�, 
                    % other constellations range between 180 and 360�
in = degtorad(72);  % inclination of the planes
s = 21;             % satellites per plane
p = 9;              % number of planes
f=8;                % f = 1.25*p*cos(i); % parameter defined graphically
eo=degtorad(20);    % elevation angle [rad]

%____Numerical Input_________
r = 1;    %radius of the sphere
n = 50;   %number of cells of the sphere generated
N = 50;   %number of points used to print the orbit
nn = 100; %number of points used to compute the footprint

%____Physical Constants_____
global Re, Re = 6371;  %Earth's Radius [Km]
%% 1. Preprocessing

global m, m = degreegen/180; 
            % generates Walker Delta Constellation (m=2 -> 2*180�
            % generated constellation), Semi Walker (m=1 -> 1*180�) 
            % or other 1<m<2.
                   
global a, a = (Re+h)/Re; 
            %radius of the orbit (in terms of the sphere radius) R=h*r;

%% 2. Main Process

%compute the sat RAAN and some sat. parameters (Laura)
Omega = zeros(1,p);
nu = zeros(s,p);     %true anomaly of the satellites
angle = zeros(1,p);  %angle between satellites
for i = 1:p
    Omega(i) = m*pi*(i-1)/p;   %orbits
    
    angle(i) = f*2*pi*(i-1)/(s*p); %sats 
    for k = 1:s
        nu(k,i) = 2*pi*(k-1)/s+angle(i);
    end
end

% Compute the satellite coordinates (Laura)
X = zeros(3,s,p);
for i = 1:p
    for k = 1:s
        X(:,k,i) = cartesian(a,0,in,Omega(i),0,nu(k,i));
    end
end
sCOOR=reshape(X,3,p*s)'; %put the data in a more friendly structure

%compute the radius of the footprint 
%SATELLITE_FOOTPRINT_COMPUTATION.m
d=Re*((((h+Re)/Re)^2-cos(eo)^2)^.5-sin(eo));
Bo=asin(d*cos(eo)/(h+Re));
S=2*pi*Re^2*(1-cos(Bo));
R=(S/pi)^.5;

%create the footprint of each satellite
SATfp=[];
for i=1:s*p
    %unitary vector pointing to the sat from the earth's center
    u=-sCOOR(i,:)/norm(sCOOR(i,:));
    %Generate the circle coordinates (2D)
    [satfp(:,1),satfp(:,2),satfp(:,3)]=circle3D(sCOOR(i,:),u, R/Re, nn);
    satfp=satfp+repmat(u*h/Re,nn,1); %move the circle to be tg to Earth
    
    for j=1:size(satfp,1)
        v=sCOOR(i,:)-satfp(j,:); v=v/norm(v); 
        %unitary vector sat->footprint
        
        %verify if the footprint intersect with the earth's surface
        tol=satfp(j,1)^2+satfp(j,2)^2+satfp(j,3)^2-r^2;
        while tol>1e-2
            satfp(j,:)=satfp(j,:)-0.01*v;
            tol=satfp(j,1)^2+satfp(j,2)^2+satfp(j,3)^2-r^2;
        end
    end
    %store the footprint
    SATfp=[SATfp; satfp];
end


%% 3. Plots

%3D plot
[COOR]=Orbitplot3D( in,r, X,n, SATfp, nn );
%Ground track
x=zeros(size(COOR,1),size(COOR,2)/3); y=x; z=x;
for i=1:size(COOR,2)/3
    x(:,i)=COOR(:,3*(i-1)+1);
    y(:,i)=COOR(:,3*(i-1)+2);
    z(:,i)=COOR(:,3*(i-1)+3);
end
GroundTrackPlot( x, y, z, sCOOR, SATfp, nn);