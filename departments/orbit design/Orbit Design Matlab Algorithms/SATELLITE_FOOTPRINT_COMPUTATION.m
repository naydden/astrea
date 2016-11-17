% COMPUTATION OF A SATELLITE FOOTPRINT

clc 
clear all
close all

%% PHISICAL CONSTANTS AND PARAMETERS 

Re = 6378;          %Earth's Radius [Km]


Se = 4*pi*Re^2;     %Earth's Surface [Km^2]

h = 500:100:1000;   %Satellte height [Km]
eo = 0:5:40;      %Elevation angle [º]

%% SOLVER

for i = 1:length(h)
    for j = 1:length(eo)
        
        d(i,j)= Re * ( sqrt( ( (h(i)+Re)/Re )^2 - ( cosd(eo(j)) )^2 ) - sind(eo(j)) );
        Bo(i,j) = asind( d(i,j)*cosd(eo(j)) / (h(i)+Re) );
        
        S(i,j) = 2 * pi * Re^2 * ( 1-cosd(Bo(i,j)) );
        
        %R(i,j) = sqrt( d(i,j)^2 - h(i)^2 );
        R (i,j) = sqrt ( S(i,j)/pi );
        
        %Satellite coverage expressed as a fraction of Earth's Area (%)
        Cov(i,j) = ( S(i,j) / Se )*100; 
        
        %Number of satellites needed
        nsat(i,j) = ( Se / S(i,j) );

    end
end

%% PLOTS

figure
for k = 1:length(h)
    plot ( eo , Cov(k,:) )
    hold on
end
title('Coverage vs Elevation')
xlabel('Elevation [deg]') % x-axis label
ylabel('Coverage [%]') % y-axis label
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)))


figure
for k = 1:length(h)
    plot ( eo , nsat(k,:) )
    legend ('h(j)')
    hold on
end
title('Num.Satellites vs Elevation')
xlabel('Elevation [deg]') % x-axis label
ylabel('Num.Satellites') % y-axis label
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)))       

figure
for k = 1:length(h)
    plot ( eo , 2*Bo(k,:) )
    hold on
end
title('Angle vs Elevation')
xlabel('Elevation [deg]') % x-axis label
ylabel('Angle between satellites[deg]') % y-axis label
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)))



% Calcular radio
% Buscar angle elevation