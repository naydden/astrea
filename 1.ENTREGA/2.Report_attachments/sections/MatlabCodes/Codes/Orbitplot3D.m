function [COOR] = Orbitplot3D( in,r, X,n, sCOOR, N )
%Orbitplot3D is a function that creates a 3D representation of a
%constellation and plots the satellites. It can be also configured to plot
%their footprint
%In:
%   in: inclinaton angle [rad]
%   r: sphere radius [-]
%   X: Coordinates of the sats (3,spp,p)
%   n: number of points used to draw thhe sphere
%   ~~o~~
%   COOR: coordinates of the footprints
%   N: number of points used to compute the footprint
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p=size(X,3);
s=size(X,2);
global a m;

% Earth
%Create the earth instance as an sphere
[x,y,z]=sphere(n); %Generate the sphere
surf(r*x,r*y,r*z,'FaceColor','interp','EdgeColor','none'); %plot the sphere
colormap(linspace(0,1,n)'*[0 0.6 1]); %custom colormap (blue tones);
axis equal;
box off;
hold on;
rotate3d on;

% Orbits
%Generate the circle coordinates (2D)
coor=zeros(N+1,3);
for i=1:N
    coor(i,:)=[a*r*cos(2*i*pi/N), a*r*sin(2*i*pi/N), 0];
end
coor(N+1,:)=coor(1,:);

%Rotation Matrices
Ry=[ cos(in) 0 sin(in)
       0     1   0   
    -sin(in) 0 cos(in)];

Rz= @(Omega) [-sin(Omega) cos(Omega) 0
              cos(Omega)  sin(Omega) 0
                  0           0      1]; % he modificat aquesta matriu
              
coor=coor*Ry; %Rotate the circle in the y axis (inclination)

%rotate the plane to get the other planes
COOR=zeros(N+1,3*p);
for i=1:p
    COOR(:,(1+(i-1)*3):(3+(i-1)*3))=coor*Rz(m*(i-1)*pi/p);
end
for i=1:p
    plot3(COOR(:,1+(i-1)*3),COOR(:,2+(i-1)*3),COOR(:,3+(i-1)*3),...
        'linewidth',1.5);
end

% Satellites (Laura)
for i = 1:p
    for k = 1:s
      scatter3(X(1,k,i),X(2,k,i),X(3,k,i),25,'square','k','filled');
      hold on;
    end
end

%Plot also the footprint if the number of input variables is bigger than 4
if nargin > 4
    for i=1:p*s
        b=1+N*(i-1); c=N+N*(i-1);
        fill3(sCOOR(b:c,1),sCOOR(b:c,2),sCOOR(b:c,3),[1 0 0],...
            'FaceAlpha','.3','LineStyle','none');
    end
end
end

