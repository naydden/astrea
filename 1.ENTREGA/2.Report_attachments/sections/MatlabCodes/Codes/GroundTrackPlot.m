function [] = GroundTrackPlot( x, y, z, coor, fp_coor, N)

%GroundTrackPlot is used to plot the ground 
%track of the constellation orbits
%given their coordinates (cartesian).It has two variants, only introducing
%the first 3 inputs will plot the orbits and if you introduce all the
%inputs it plots the satellites and their footprint.
%In:
%   x: matrix containing the x coordinate of the orbits to plot 
%   y: matrix containing the y coordinate of the orbits to plot 
%   z: matrix containing the z coordinate of the orbits to plot 
%   ~~o~~
%   coor: satellite coordinates (nx3)
%   fp_coor: satellite footprint coordinates ((n*nn)x3)
%   N: number of points used to compute the footprint
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p=size(x,2); %n of planes
s=size(coor,1)/size(x,2); %n of spp

%% Set the background image for the ground track plot
figure('name','Constellation Ground Track');
tit=['Num. of planes: ' num2str(p) ' Num. of spp: ' num2str(s)];
ground='World-satellite-map.png';
im=imread(ground);
imagesc([-180,180],[-90,90],im);
grid on;
axis([-180 180 -90 90]);
title(tit);
hold on;
%set the plot lines color
c = ['r','y','m','g','w']; C=[];
for i=1:ceil(size(x,2)/5)
    C=[C c];
end

%% Main process
for i=1:size(x,2)
    %convert the cartesian coordinates to spherical
    [lon,lat]=cart2sph(x(:,i),y(:,i),z(:,i)); 
    %get the coordinates of the orbit in spherical coord.
    lon=radtodeg(lon); %pass to degrees
    lat=radtodeg(lat);
    Lon=[lon(end); lon(1:end-1)];
    
    %Plot the sat ground track
    [~,j]=max(abs(lon-Lon));
    plot(lon(1:j-1),lat(1:j-1),C(i),'linewidth',1.2);
    plot(lon(j:end),lat(j:end),C(i),'linewidth',1.2);
end

if nargin > 3
    [slon,slat]=cart2sph(coor(:,1),coor(:,2),coor(:,3)); 
    %get the coordinates of the sats in spherical coord.
    slon=radtodeg(slon); %pass to degrees
    slat=radtodeg(slat);
    
    %plot the sats
    scatter(slon,slat,30,'square','w','filled');
    
    %get the coordinates of the sats footprint in spherical coord.
    [fplon,fplat]=cart2sph(fp_coor(:,1),fp_coor(:,2),fp_coor(:,3)); 
    fplon=radtodeg(fplon); %pass to degrees
    fplat=radtodeg(fplat);
    
    %Plot their footprint
    for i=1:p*s
        b=1+N*(i-1); c=N+N*(i-1);
        fpLon=[fplon(c); fplon(b:c-1)];
        
        %get the maximum values (plotting purposes)
        [sorted, I]=sort(abs(fplon(b:c)-fpLon));
        M=sorted(end-1:end);
        j=I(end-1:end); j(2)=j(2)-1;
        
        if M(1) > 300
            if j(1) < j(2)
                auxlon=fplon(b+j(1):b+j(2)-1); %negative lon. part 
                auxlat=fplat(b+j(1):b+j(2)-1);

                auxlon2=[fplon(b+j(2):c); fplon(b:b+j(1)-2)]; %+ lon. part
                auxlat2=[fplat(b+j(2):c); fplat(b:b+j(1)-2)];
            else
                auxlon2=fplon(b+j(2):b+j(1)-2); %negative lon. part
                auxlat2=fplat(b+j(2):b+j(1)-2);

                auxlon=[fplon(b+j(1):c); fplon(b:b+j(2)-1)]; 
                %positive lon. part
                auxlat=[fplat(b+j(1):c); fplat(b:b+j(2)-1)];
            end
            
         fill(auxlon,auxlat,[1 0 0],'FaceAlpha','.2','LineStyle','none');
         fill(auxlon2,auxlat2,[1 0 0],'FaceAlpha','.2','LineStyle','none');
        else
            fill(fplon(b:c),fplat(b:c),[1 0 0],...
                'FaceAlpha','.3','LineStyle','none');
        
        %%%%%%%%%%%%
        pause(.05); %
        %%%%%%%%%%%%
        end
        
        
    end
end
end