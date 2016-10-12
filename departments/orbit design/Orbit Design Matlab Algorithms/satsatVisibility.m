function [phi, d] = satsatVisibility(hatm, h)
%This matlab routine computes the maximum distance between satellites only 
%considering that they must see each other. We neglect the attenuation effects
%of the atmosphere (the only restriction is to avoid the first hatm km)
%In: 
%  hatm = estimated height (atmosphere) for fully-absortion of the signal  [km]
%  h = satellite height [km]
%Out:
%  phi = angle between the satellites in the same plane [rad]
%    d = maximum distance between satellites [km]

R = 6371.001; %earth Radius [km]

%compute the angle between sattellites
phi = 2*acos((R+hatm)/(R+h));
d = (2*(R+h)^2-2*(R+h)^2*cos(phi))^.5;
end
