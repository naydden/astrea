% This function converts kepler orbit elements to cartesian coordinates

% a: semimajor axis
% e: eccentricity
% i: inclination [rad]
% Omega: longitude of the ascending node [rad]
% w: Argument of periapsis [rad]
% nu: True anomaly [rad]

function [y] = cartesian2(a,e,i,Omega,w,nu)
global mu;
% Position in cylindrical coordinates
r = a*(1-e^2)/(1+e*cos(nu));

% Position components
X = [r*(cos(Omega)*cos(w+nu)-sin(Omega)*sin(w+nu)*cos(i));
    r*(sin(Omega)*cos(w+nu)+cos(Omega)*sin(w+nu)*cos(i));
    r*sin(i)*sin(w+nu)];

A = sin(Omega)*cos(i)*(e*cos(w)+cos(w+nu));
v(1)=-sqrt(mu/( a*(1-e^2)))*(cos(Omega)*(e*sin(w)+sin(w+nu))+A);

B = cos(Omega)*cos(i)*(e*cos(w)+cos(w+nu));
v(2)=-sqrt(mu/( a*(1-e^2)))*(sin(Omega)*(e*sin(w)+sin(w+nu))-B);

v(3)=sqrt(mu/( a*(1-e^2)))*(e*cos(w)+cos(w+nu))*sin(i);

y=[X;v'];
end