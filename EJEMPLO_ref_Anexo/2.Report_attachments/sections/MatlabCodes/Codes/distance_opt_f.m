function [anglesp] = distance_opt_f (i,a,p,s,typeWD,F)
%% Data input

% i: Inclination angle [deg]
% a: Height of the orbit [km]
% p: Number of planes
% s: Number of satelites per plane

%% Walker Delta type of geometry
WD=[180 225 360];
degreegenerate = WD(typeWD); 
% Walker Delta 360º, SemiWalker Delta 180º, 
% other constellations range between 180 and 360�

m = degreegenerate/180; 
% generates Walker Delta Constellation 
% (m=2 -> 2*180� generated constellation), 
% Semi Walker (m=1 -> 1*180�) or other 1<m<2.

% Phasing between adjacent planes
%f = F*p*cosd(i); % parameter defined graphically
 f = F; % Whereas F is a number 0 < f < (p-1)
        % According to Astrodynamics notes
        
REarth = 6371; % [km]
h = (REarth+a)/REarth; %radius of the orbit 
%(in terms of the sphere radius) R=h*r;

%% Distribution of Coordinates of the satellites
Omega = zeros(1,p);
nu = zeros(s,p);
angle = zeros(1,p);
X = zeros(3,s,p);

for j = 1:p
    Omega(j) = m*pi*(j-1)/p;
    angle(j) = f*2*pi*(j-1)/(s*p); %Phasing due to f between planes
    for k = 1:s
        % True anomalies of the s satellites
        nu(k,j) = 2*pi*(k-1)/s+angle(j);
        X(:,k,j) = cartesian(h,0,i,Omega(j),0,nu(k,j));
        % h[adim] --> X[adim]
    end
end

%% Fast comprovation
% Does the adjacent phasing exceed pi 
% THIS MEANS
% The minimum might be between 


%% Angle between different satellites

c=zeros(s,p); %General aproach. Where is the minimum angle?)

for j=2:p
  for i=1:s-1
        
    % We need to assess the angles between the last and first plane
    % Specially in Semi-Walker configurations. That's why we add this
    % auxiliar variables.
        
      p1=j;
      if j==1
          p2=p;
      else
          p2=j-1;
      end
        
    % Angles a2 and a3 could be used to increase de reliability of the 
    % constellation, allowing greater possibilities of communication 
    % between planes. (not only between (i,p1) and (i,p2)
        
    % Angle between satellite (i,p1) and satellite (i,p2)
    a1=acosd((X(:,i,p1)'*X(:,i,p2))/(norm(X(:,i,p1))*norm(X(:,i,p2))));
        
    % Angle between satellite (i,j) and satellite (i+1,j-1)
    a2=acosd((X(:,i,p1)'*X(:,i+1,p2))/(norm(X(:,i,p1))*norm(X(:,i+1,p2))));

    % Angle between satellite (i,j) and satellite (i-1,j-1)
    a3=acosd((X(:,i,p1)'*X(:,i-1,p2))/(norm(X(:,i,p1))*norm(X(:,i-1,p2))));

    % Among the computed angles, choose the minimum
    [angle(i,j),c(i,j)]=min([a1 a2]);
        
    % Why minimum? Then we know that at least it will be able to 
    % cover the angle with one of the two adjacent close nodes.
        
  end
end

anglesp = max(max(angle));             %Maximum angle between satellites of different planes
