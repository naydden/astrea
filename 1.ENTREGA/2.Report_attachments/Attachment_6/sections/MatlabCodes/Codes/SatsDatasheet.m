%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% SATELLITES DATASHHET

% Writing the constellation caracteristics

a=6378.01e3;
u=3.986012e14;
h=(a+542e3)/a; 
ap=h*a;
J2=1.0826e-3;% J2=0;

ws=sqrt(u/ap^3); % Satellite mean motion
i=72*pi/180;

p=9;
spp=21;
N=p*spp;
%f=floor(0.25*spp);
f=0;
fprintf('\ni=%gdeg, p=%g, spp=%g, N=%g\n',i*180/pi,p,spp,N);


%Omega = 0:30:240; Omega=Omega*pi/180;
Omega = 0:225/(p-1):225; %Omega=Omega*pi/180;
%Omega=0:360/p:360-360/p; Omega=Omega*pi/180;

P=2*pi/ws/60;

%% Excel IMPORT

D=zeros(N,9);
for sat=1:N
    
    Om=ceil(sat/spp); ID_p=sat-(Om-1)*spp;
    ph=(ID_p-1)*360/spp;
   %[ID Plane h P i e Omega phase per]
    D(sat,1:9)=[sat Om (ap-a)/1000 P i*180/pi 0 Omega(Om) ph 0];
    
end
NameCells=['B2:I' num2str(N)];
xlswrite('SatsDatasheet.xls',D,NameCells)

%% Latex IMPORT

for sat=1:N
    
    Om=ceil(sat/spp); ID_p=sat-(Om-1)*spp;
    ph=(ID_p-1)*360/spp;
   %[ID Plane h P i e Omega phase per]
    D(sat,1:9)=[sat Om (ap-a)/1000 P i*180/pi 0 Omega(Om) ph 0];
    
end

[rows,cols]=size(D);
file=fopen('Autotable.txt','w');

for i=1:rows
    Text=['AstreaSAT ' num2str(i) ' $ '];
    for j=2:cols
        if j~=cols
            Text=[Text num2str(D(i,j-1)) ' $ '];
        else
            Text=[Text num2str(D(i,j-1)) ' \\ '];
        end
    end
    fprintf(file,[Text '\n']);
end

fclose(file);