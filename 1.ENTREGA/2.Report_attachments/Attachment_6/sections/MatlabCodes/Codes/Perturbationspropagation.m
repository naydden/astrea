%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% ORBIT PERTURBATIONS - Influential phenomena computation

%PROBLEM:
% Given: - Initial orbital parameters
% Compute the final orbtial parameters after each orbit
%clear; clc;

%% 1.A Input data physical constants
% PHISICAL CONSTANTS AND PARAMETERS 
%Physical data
RE=6.378e6; %Earth Radius [m]
u=3.986e14; %GM Earth
J2=1082.6e-5; %J2 coefficient;

%% 1. Input data Orbit parameters
% SPACECRAFT DATA
m=4;        % Satellite mass [kg]
Cd=2.2;     % Drag Coefficient
A=0.1;  % Satellite surface [m^2] --> 3U pointing to Earth
Bc=m/(Cd*A);% Ballistic coefficient of the satellite; 

% ORBIT DATA
H0=542e3;     % Altura inicial [m]
e0=0.0001;     % Excentricity
i0=72*pi/180;       % Inclination [deg]
a0=RE+H0;

% INITIAL PARAMETERS
P=2*pi*(a0^3/u)^.5; % Orbit Period [s]
v=sqrt(u/a0);       %satellite velocity [m/s]
n=(2*pi/P)*86400;  %number of revolutions/day


% SIMULATION PARAMETERS
N=10000000;

%% 2. Perturbations propagation
%Outputs de la funcio Perturbation:
%   Pert: matrix (4x5) with the perturbations causes at the rows (J2,Drag,
%   Moon,Sun) and the orbital elements in the columns (a,e,i,omega,Omega)
%     a  e  i  w  Omega  
%   [                   ]J2
%   [                   ]Drag
%   [                   ]Moon
%   [                   ]Sun as 3rd Body
 

    
t=zeros(1,N); t(1)=0;
a=t; a(1)=a0;
e=t; e(1)=e0;
i=t; i(1)=i0;

w_sun=zeros(1,N-1); w_moon=w_sun; w_J2=w_sun;
Om_sun=w_sun; Om_moon=w_sun; Om_J2=w_sun;

w0=0; Omega0=0;
w=t; w(1)=w0;
Omega=t; Omega(1)=Omega0;

tic
for n=2:N
    
    a0=a(n-1);
    e0=e(n-1);
    i0=i(n-1);
    w0=w(n-1);
    Omega0=Omega(n-1);
    
    P=2*pi*(a0^3/u)^.5; % Orbit Period [s]
    
    [Perti rho]=Perturbation(a0,e0,i0,Bc);
    
    a(n)=a0+sum(Perti(:,1));
    e(n)=e0+sum(Perti(:,2));
    i(n)=i0+sum(Perti(:,3));
    w(n)=w0+sum(Perti(:,4)); 
    Omega(n)=Omega0+sum(Perti(:,5)); 
    
    w_sun(n-1)=Perti(4,4);
    w_moon(n-1)=Perti(3,4);
    w_J2(n-1)=Perti(1,4);
    
    Om_sun(n-1)=Perti(4,5);
    Om_moon(n-1)=Perti(3,5);
    Om_J2(n-1)=Perti(1,5);    
    
    % We don't want angles bigger than 360!
    w(n)=w(n)-360*floor(w(n)/360);
    Omega(n)=Omega(n)-360*floor(Omega(n)/360);
    
    if Omega(n)>2*pi
        Omega(n)=Omega(n)-2*pi;
    end
    
    t(n)=t(n-1)+P;
    
    if a(n)<(RE+180e3)
     fprintf('Your satellite succesfully burned in the atmosphere.\n\n')
      break
    end
    
end

%% Post processing and results plotting

time=toc;
fprintf('Time to do %g iterations = %g s.\n ',n,time)

figure(1)
plot(t(1:n)/(3600*24),(a(1:n)-RE)/1000)
ylim([100 500])
grid on
ylabel('Orbit height [km]')
xlabel('Time [days]')
tfin=floor(t(n)/(3600*24));
days = num2str(floor(tfin/365));
titulaso=['Orbit decay in ' num2str(tfin) ' days = ' days 'years']; 
title(titulaso)
 

figure(2)
subplot(1,2,1)
plot(t(1:n)/(3600*24),w(1:n))
grid on
ylabel('Perigee Argument [deg]')
xlabel('Time [days]')
title('Perigee Argument deviation in 100 days')
axis([0 100 0 360])

subplot(1,2,2)
plot(t(1:n)/(3600*24),Omega(1:n))
grid on
ylabel('Ascenent Node Argument [deg]')
xlabel('Time [days]')
tfin=t(n)/(3600*24);
title('Ascendent Node deviation in 100 days')
axis([0 100 0 360])

%% Modulus Analysis

figure(3)

subplot(1,2,1)
semilogy(t(2:n+1)/(3600*24),abs(w_sun(1:n)),'b','LineWidth',2);
hold on
semilogy(t(2:n+1)/(3600*24),abs(w_moon(1:n)),'r','LineWidth',2)
semilogy(t(2:n+1)/(3600*24),abs(w_J2(1:n)),'k','LineWidth',2)
legend('Sun','Moon','J2')
xlabel('Time (days)');
ylabel('Modulus of Perturbation []');
title('Effects on Argument of the Perigee');
grid on
axis([0 t(n)/(3600*24) 1e-8 1])

subplot(1,2,2)
semilogy(t(2:n+1)/(3600*24),abs(Om_sun(1:n)),'b','LineWidth',2);
hold on
semilogy(t(2:n+1)/(3600*24),abs(Om_moon(1:n)),'r','LineWidth',2)
semilogy(t(2:n+1)/(3600*24),abs(Om_J2(1:n)),'k','LineWidth',2)
legend('Sun','Moon','J2')
xlabel('Time (days)');
ylabel('Modulus of Perturbation []');
title('Effects on Argument of Ascendent Node');
axis([0 t(n)/(3600*24) 1e-8 1])
grid on

print -depsc ModulusAngulars