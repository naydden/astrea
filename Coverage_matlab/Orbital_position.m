function [X] = Orbital_position(N_sat,N_planes,h,I,mu,t)
%N_sat: number of sats per plane
%N_planes: number of planes
%h: Orbit high in km
%I: inclination in degrees
%mu: Interplane desphase in degrees
%t: array of time in minutes

I=I*pi/180;   %I in rad
mu=mu*pi/180; %mu in rad
Mu=0:mu:(N_planes-1)*mu; %array of absolut desphase of the planes

L_1o=[1 0 0; 0 cos(I) -sin(I); 0 sin(I) cos(I)]; %1st rotation matrix

Re=6371;    %earth radius km
h=h/Re;     %high in earth radius
GM=6.474e-11*5.972e24*3600/Re^3; %Gravitational constant*earth mass [(Re^3)/min^2]
Ro=1+h;     %orbital radius in earth radius

w=sqrt(GM/(Ro^3));  %sat's angular speed [rad/min]
phi=0:2*pi/N_sat:(N_sat-1)*2*pi/N_sat; %absolut desphase of each sat of a plane respect the first 
plane=1:1:N_planes;
f=0;
phi_0=f*2*pi*(plane-1)/(N_planes*N_sat);
l=length(t);

for i=1:N_planes
    L_e1=[cos(Mu(i)) -sin(Mu(i)) 0; sin(Mu(i)) cos(Mu(i)) 0; 0 0 1]; %2nd rotation matrix
    L_eo=L_e1*L_1o;   %Global rotation matrix orbit coordinates to global coordinates
    for j=1:N_sat
        for k=1:l
            X_loc=[Ro*cos(w*t(k)+phi(j)+phi_0(i));0;-Ro*sin(w*t(k)+phi(j)+phi_0(i))]; %orbit coordinates
            X(:,(i-1)*N_sat+j,k)=L_eo*X_loc;    %global coordinates
        end
    end        
end
end

