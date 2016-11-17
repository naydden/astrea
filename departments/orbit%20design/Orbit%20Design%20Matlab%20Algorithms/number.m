% Program that calculates the number of satellites necessary to get global
% coverage in a LEO for different constellations

RE = 6371e3; % Earth radius [m]
REquator = 6378137; % [m]
RPolar = 6356752.3; % [m]

% Altitude [m]
h = 300e3:100e3:1000e3;
hatm = 200e3;

% Elevation angle [rad]
fepsilon = 20*pi/180;
epsilon = linspace(0,fepsilon);
epsilonmin = epsilon;

%% Polar

t = zeros(length(h),length(epsilon));

for i = 1:length(h)
    
    rs = RE+h(i); % radius of the orbit
    
    % Angular radius of the Earth [rad]
    rho = asin(RE/(RE+h(i)));
    
    % Maximum distance between satellites
    [maxtheta, ~] = satsatVisibility(hatm/1000,h(i)/1000);
    
    for j = 1:length(epsilon)
        
        % Maximum nadir angle
        etamax = asin(sin(rho)*cos(epsilonmin(j)));
        
        % Maximum Earth central angle
        lambdamax = pi/2-epsilonmin(j)-etamax;
        
        arcEq = 2*REquator*lambdamax;
        p = ceil(2*pi*REquator/arcEq); % number of planes
        
        s = ceil(2*pi/maxtheta); % satellites per plane
        
        t(i,j) = p*s;
    end
end

figure(1);
for i = 1:length(h)
    plot(epsilon.*180/pi,t(i,:));
    hold on
end
legend('300','400','500','600','700','800','900','1000');
xlabel('Elevation angle (º)');
ylabel('Number of satellites');
title('Polar constellation');

%% Walker delta

% Walker delta parameters
% t = total number of satellites
% p = orbital planes
% s = satellites per plane
% f = relative spacing between planes

