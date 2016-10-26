%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MINIMUM INCLINATION - To provide Full coverage %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PROBLEM:
% Given: - Minimum Elevation over the horizon to contact the satellite
%        - Height of the orbit
% Compute the minimum elevation to ensure visual contact
clear; clc;

%% Key variable
min_elevation_pole=5; %Initially set by @lfore as 32;

%% Previous Calculation: Visibility(latitude)
% Extract from ElevationAngle.m
% S/O to @lfore! Good job Team

R = 6371; 
N = 180;
deg = 180;
x = linspace(0,90,N);

%%%% Minimum angle of elevation due to atmospheric conditions at given
%%%% latitudes. 0 km represents the poles whereas 6371 km represents the
%%%% Earth's equator taking as reference point the south pole.

elvlat = zeros (1,deg);
for i = 1:deg/5    
    elvlat(i) = min_elevation_pole-(min_elevation_pole-15)/(2*(90-75))*i;
end
for i = deg/5+1:4*deg/10    
    elvlat(i) = elvlat(deg/5);
end
for i = 4*deg/10+1:6*deg/10    
    elvlat(i) = elvlat(2*deg/5)+3/15*(i-4*deg/10);
end
for i = 6*deg/10+1:8*deg/10
    elvlat(i) = elvlat(6*deg/10)-3/15*(i-6*deg/10);
end
for i = 8*deg/10+1:10*deg/10    
    elvlat(i) = elvlat(8*deg/10)-1/15*(i-8*deg/10);
end
elvlat=fliplr(elvlat); %Sets the first value to the value in the equator
                        %The first value matches 0 latitude

%% Minimum i computation

h = 400:100:1000;   %Satellte height [Km]
L = 0:5:90;         %Latitude angle [º]
e = (interp1q(x',elvlat',L'))'; %Minimum elevations interpolation [º]

figure(1)
subplot(1,2,1)
plot(x,elvlat)
subplot(1,2,2)
title('Interpolation Verification')
plot(L,e)

inc=zeros(length(h),length(L)); %Inclinations preallocation

for i = 1:length(h)
    for j = 1:length(L)
               
        A=cosd(e(j))/(1+h(i)/R);
        theta=acosd(A)-e(j); %Earth Central Angle
        inc(i,j)=L(j)-theta;
        if inc(i,j)<0
            inc(i,j)=0;
        end
    end
end

%% 3. Results Plotting
figure(2)
colors=['r','b','y','g','m','k','c'];
for k = 1:length(h)
    hold on
    plot ( L , inc(k,:),colors(k))
end
grid on
title('Minimum inclination required to fulfill global coverage')
xlabel('Latitude [º]')
ylabel('Inclination [º]')
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)),num2str(h(7))) 

