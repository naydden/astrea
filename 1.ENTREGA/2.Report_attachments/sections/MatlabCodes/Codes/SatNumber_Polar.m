%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% POLAR ORBITS Number of satellites computation

%PROBLEM:
% Given: - Elevation over the horizon to contact the satellite
%        - Height of the orbit
% Compute the final number of satellites
% Through the "Streets of coverage" method.

clear; clc;

%% 1. Input data
% PHISICAL CONSTANTS AND PARAMETERS 

Re = 6378;          %Earth's Radius [Km]
Se = 4*pi*Re^2;     %Earth's Surface [Km^2]
h = 500:20:600;   %Satellte height [Km]
eo = 20;        %Elevation angle [�]

%% 2. Number of satellites computation

d=zeros(length(h),length(eo));
nsats_opt=d;

for i = 1:length(h)
    for j = 1:length(eo)
        
        d(i,j)= Re*(sqrt(((h(i)+Re)/Re)^2-(cosd(eo(j)))^2)-sind(eo(j)));
        Bo = asind( d(i,j)*cosd(eo(j)) / (h(i)+Re) );
        Sup = 2 * pi * Re^2 * ( 1-cosd(Bo) );
  
        Nspp_min=ceil(360/(2*Bo));        %Number of satellites per plain
        Nsatspp=Nspp_min:1:Nspp_min+50;   %Array of sats/plain to optimize
        
        S=360./Nsatspp;
        
        Lstreet=acosd(cosd(Bo)./cosd(S/2));
        
        DmaxS=Lstreet+Bo;
        DmaxO=2*Lstreet;

        Nplains=ceil(1+((180-DmaxO)./DmaxS));
        Nsats=Nplains.*Nsatspp;
        [row,col]=find(min(Nplains));
        Nplains_opt(i,j)=Nplains(row,col);
        nsats_opt(i,j)=Nsats(row,col);
        
        %Number of satellites needed
        nsat(i,j) = ( Se / Sup );
        
    end
end

%% 3. Results Plotting

figure(2)
colors=['r','b','y','g','m','k','c'];
for k = 1:length(eo)
    subplot(1,2,1)
    hold on
    semilogy ( h , nsats_opt(:,k),colors(k))
    subplot(1,2,2)
    hold on
    semilogy ( h , Nplains_opt(:,k),colors(k))
end

subplot(1,2,1)
title('Num.Sats vs Height')
xlabel('Height [km]') % x-axis label
ylabel('Num.Satellites') % y-axis label
%legend(num2str(eo(1)),num2str(eo(2)),num2str(eo(3)),...
%    num2str(eo(4)),num2str(eo(5)),num2str(eo(6)),num2str(eo(7)))
axis([400 900 0 1000])
grid on

subplot(1,2,2)
title('Num.Planes vs Height')
xlabel('Height [km]') % x-axis label
ylabel('Num of Orbital Planes') % y-axis label
% legend(num2str(eo(1)),num2str(eo(2)),num2str(eo(3)),...
%     num2str(eo(4)),num2str(eo(5)),num2str(eo(6)),num2str(eo(7)))   
axis([400 900 0 30])
grid on

%% Single detailed analysis

figure(1)
plot(h,nsats_opt)
title('Num.Sats vs Height - Streets of coverage Method')
xlabel('Height [km]') % x-axis label
ylabel('Num.Satellites') % y-axis label
grid on

%% Multianalysis variation with height

figure(3)
colors=['r','b','y','g','m','k'];
for k = 1:length(h)
    subplot(1,2,1)
    hold on
    plot ( eo , nsats_opt(k,:),colors(k))
    subplot(1,2,2)
    hold on
    plot ( eo , Nplains_opt(k,:),colors(k))
end

subplot(1,2,1)
title('Num.Sats vs Elevation - Streets of coverage Method')
xlabel('Elevation [deg]') % x-axis label
ylabel('Num.Satellites') % y-axis label
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)))   
axis([0 40 0 700])
grid on

subplot(1,2,2)
title('Num.Planes vs Elevation - Surfaces Method')
xlabel('Elevation [deg]') % x-axis label
ylabel('Num.Orbital Planes') % y-axis label
legend(num2str(h(1)),num2str(h(2)),num2str(h(3)),...
    num2str(h(4)),num2str(h(5)),num2str(h(6)))  
axis([0 40 0 30])
grid on