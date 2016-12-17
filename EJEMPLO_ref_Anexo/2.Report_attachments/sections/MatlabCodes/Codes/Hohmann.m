function [deltaV1,deltaV2,mp,t] = Hohmann(hinicial,hfinal,m,Isp)

% Hohmann transfer orbit between two circular orbits
% - hinicial: height of the first orbit [m]
% - hfinal: height of the second orbit [m]
% - m: mass of the satellite (total) [kg]
% - Isp: Specific impulse [s]
% - mp: fuel mass [kg]
% - t: time needed to do the maneuver [s]

mu = 3.986004418e14; % Standard gravitational parameter (Earth)
REarth = 6.371e6; % [m]
g0 = 9.81; % Earth's gravity [m/s^2]

% First orbit
r1 = REarth+hinicial; % [m]
v1 = sqrt(mu/r1); % [m/s]

% Second orbit
r2 = REarth+hfinal; % [m]
v2 = sqrt(mu/r2); % [m/s]

% Transfer orbit (ellipse)
vp = sqrt(2*mu*r2/(r1*(r1+r2))); % [m/s]
va = sqrt(2*mu*r1/(r2*(r1+r2))); % [m/s]
a = (r1+r2)/2;
T = sqrt(2*pi^2*a^3/mu);

% deltaV 1 -> p(transfer orbit)
deltaV1 = vp-v1;
% deltaV a(transfer orbit) -> 2
deltaV2 = v2-va;
% Total
deltaV = deltaV1+deltaV2;

% Fuel mass
mp = m*(1-exp(-deltaV/(g0*Isp)));

% Time
t = T/2;

end