%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% ORBIT SMART WALKER-DELTA

% This routine compute the WD possible configurations that give global
% coverage by having the distance between planes matches the planes
% rotation

p=[12 15 18];
S=360./p; S=S*pi/180;

re=6378.01e3;
u=3.986012e14;
J2=1.0826e-3;

we=2*pi/(24*3600);

syms a;
inc=75; inc=inc*pi/180;
a_sol=zeros(1,length(inc));

for i=1:length(inc)
    dOmega=-1.5*J2*(re/a)^2*sqrt(u/a^3)*cos(inc(i));
    P0=2*pi*sqrt(a^3/u);
    Pn=P0*(1-1.5*J2*(re/a)^2-0.75*J2*(4-5*sin(inc(i))^2)*(re/a)^2);
    eq=S(2)-Pn*(we-dOmega);
    sol=vpasolve(eq,re*1.15);
    if imag(sol)==0
        a_sol(i)=sol;
    end
end

plot(inc*180/pi,(a_sol-re)/1000)
h=(a_sol-re)/1000

%% Satellites computation

emin=20*pi/180;
Lstreet=5*pi/180;

rho=re/(re+h*1000);
etha=asin(sin(rho)*cos(emin));
Lmax=pi/2-emin-etha;
S=2*acos(cos(Lmax)/cos(Lstreet));
N=ceil(2*pi/S)