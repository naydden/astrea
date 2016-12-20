%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% THRUST 

% This function computes the deltaV and propellant mass necessary to
% maintain an orbit between to heights

%function [mp,deltaV,tHoh] = thrust(hmax,hmin,ms,Isp,Thr)

% This function computes the deltaV and propellant mass necessary to
% maintain an orbit between to heights

% Input variables:
% - hmax: maximum height [m]
% - hmin: minimum height [m]
% - ms: dry mass of the spacecraft [kg]
% - Isp: specific impulse of the spacecraft [s]
% - Thr: thrust of the spacecraft [N]

% Output variables:
% - mp: array of propellant mass necessary for every Hohmann transfer [kg]
% - deltaV: array of deltaVs necessary for every Hohmann transfer [m/s]
% - tHoh: array of time necessary to do a Hohmann transfer [s]

% Proposed values:
 hmax = 550e3;
 hmin = hmax-8;
 ms = 3.95;       % Dry mass [kg]
 Isp = 2150;
 Thr = 100e-6;

%% Data

% PHYSICAL CONSTANTS AND PARAMETERS 
RE = 6.378e6; %Earth Radius [m]
mu = 3.986e14; %GM Earth
g0 = 9.81;

% SPACECRAFT DATA
m = 4;        % Satellite mass [kg]
Cd = 2.2;     % Drag Coefficient
A = 0.1*0.3;  % Satellite surface [m^2] --> 3U pointing to Earth

% ORBIT DATA
H0 = hmax;     % Altura inicial [m]
E0 = 0.01;     % Excentricity
I0 = 80*pi/180;       % Inclination [deg]
A0 = RE+H0;

% SIMULATION PARAMETERS
N=100000;
M = 10000;

%% 2. Perturbations propagation

% Creation of the matrices
temp=zeros(M,N);
temp(1,1)=0;
H = zeros(M,N); 
H(1,1) = H0;
% deltaV = zeros(2,N);

% Asign initial values
W0=0; OMEGA0=0;
a=zeros(1,N); a(1)=A0;
e=a; e(1)=E0;
i=a; i(1)=I0;
t=a; t(1)=0;
w=a; w(1)=W0;
Omega=a; Omega(1)=OMEGA0;

for j = 1:M
    
    n = 1;
    while a(n)>=(RE+hmin)
        
        n = n+1;
        
        a0=a(n-1);
        e0=e(n-1);
        i0=i(n-1);
        w0=w(n-1);
        Omega0=Omega(n-1);
        
        P=2*pi*(a0^3/mu)^.5; % Orbit Period [s]
        
        Bc = m/(Cd*A);% Ballistic coefficient of the satellite
        Perti=Perturbation(a0,e0,i0,Bc);
        
        a(n)=a0+sum(Perti(:,1));
        e(n)=e0+sum(Perti(:,2));
        i(n)=i0+sum(Perti(:,3));
        w(n)=w0+sum(Perti(:,4));
        Omega(n)=Omega0+sum(Perti(:,5));
        
        % We don't want angles bigger than 360!
        w(n)=w(n)-360*floor(w(n)/360);
        Omega(n)=Omega(n)-360*floor(Omega(n)/360);
        
        if Omega(n)>2*pi
            Omega(n)=Omega(n)-2*pi;
        end
        
        t(n)=t(n-1)+P;
        
    end
    
    % Hohmann
    [deltaV1,deltaV2,mpit,tHohit] = Hohmann(a(n)-RE,hmax,m,Isp);
    m = m-mpit;
    
    if m<=ms
        break
    end
    
    deltaV(1,j) = deltaV1; % first row of the column -> deltaV1
    deltaV(2,j) = deltaV2; % second row of the column -> deltaV2
    mp(j) = mpit;
    tHoh(j) = tHohit;
    
    a(1)=A0;
    e(1)=E0;
    i(1)=I0;
    w(1)=W0;
    Omega(1)=OMEGA0;
    altura = a-RE;
    
    H(j,:) = altura;
    temp(j,:) = t;
    
end

%% Post process

mfr = Thr/(g0*Isp);

for i = 1:N
    if H(1,i)<=0
        break
    end
end
for j = 1:M
    if H(j,1)<=0
        break
    end
end
H = H(1:(j-1),1:(i-1));

for i = 2:N
    if temp(1,i)<=0
        break
    end
end
for j = 2:M
    if temp(j,2)<=0
        break
    end
end
temp = temp(1:(j-1),1:(i-1));
for j = 2:size(temp,1)
    temp(j,:) = temp(j,:)+tHoh(j-1)+temp(j-1,size(temp,2));
end

% Convert matrix to vector
H = reshape(H.',1,size(H,1)*size(H,2));
temp = reshape(temp.',1,size(temp,1)*size(temp,2));

figure(1)
plot(temp/(3600*24),H/1000)
grid on
ylabel('Orbit height [km]')
xlabel('Time [days]')
tfin=floor(temp(length(temp))/(3600*24));
days=num2str(floor(tfin/365));
titulaso=['Orbit decay in ' num2str(tfin) ' days = ' days ' years']; 
title(titulaso)

%% Final calculations

% % Perimeter of the ellipse
% r1 = RE+hmax; % semimajor axis
% r2 = RE+hmin; % semiminor axis
% h = (r1-r2)^2/(r1+r2)^2;
% % Ramanujan approximation
% C = pi*(r1+r2)*(1+3*h/(10+sqrt(4-3*h)));
% 
% % Check if deltaV is possible for the given thruster
% OK = true;
% for j = 1:length(mp)
%     for i = 1:2
%         mfrnecessary = deltaV(i,j)*mp(j)/(C/2);
%         if mfrnecessary>mfr
%             OK = false;
%         end
%     end
% end
% 
% if OK ==true
%     fprintf('The trajectory is possible\n\n');
% else
%     fprintf('Trajectory not possible\n\n');
% end
if length(mp)<M
    fprintf('There is still propellant left ._.\n\n')
end