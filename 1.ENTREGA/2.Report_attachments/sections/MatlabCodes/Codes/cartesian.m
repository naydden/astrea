% This function converts kepler orbit elements to cartesian coordinates

% a: semimajor axis
% e: eccentricity
% i: inclination [rad]
% Omega: longitude of the ascending node [rad]
% w: Argument of periapsis [rad]
% nu: True anomaly [rad]

function [X] = cartesian(a,e,i,Omega,w,nu)

% Position in cylindrical coordinates
r = a*(1-e^2)/(1+e*cos(nu));

% Position components
X = [r*(cos(Omega)*cos(w+nu)-sin(Omega)*sin(w+nu)*cos(i));
    r*(sin(Omega)*cos(w+nu)+cos(Omega)*sin(w+nu)*cos(i));
    r*sin(i)*sin(w+nu)];

end