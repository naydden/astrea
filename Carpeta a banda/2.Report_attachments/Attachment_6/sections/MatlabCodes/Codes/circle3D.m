function [x,y,z]=circle3D(center,normal,radius,n)
%This matlab function generates the points of a circle normal to a 3D
%vector
%In:
%   center: center of the circle (x, y, z)
%   normal: vector normal to the circle
%   radius: radius of the circle
%   n: number of points (small number can make the circle not smooth)
%Out:
%   x, y z: vectors cintaining the coord.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

theta=linspace(0,2*pi,n);
v=null(normal);
points=repmat(center',1,size(theta,2))+...
    radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
x=points(1,:); y=points(2,:); z=points(3,:);
end
