function [X] = Orbital_position(N_sat,N_planes,h,I,mu,t)
%N_sat: number of sats per plane
%N_planes: number of planes
%h: Orbit high in km
%I: inclination of the orbital plane in degrees
%mu: Angle between two adjacent orbital planes at the equatiorial plane in
%degrees
%t: array of time in minutes

I=I*pi/180;   %I in rad
mu=mu*pi/180; %mu in rad
Omega=0:mu:(N_planes-1)*mu; %longitude of the ascending node respect the inertial X axe [rad]

%for defining the orbits are defined a system of local axes X_o in which the
%y_o is the rotation axe and x_o is contained in the equatiorial plane. The
%transformation for the inertial system X_I to the local system 2 rotations are
%defined. 
%FIRST ROTATION: it defines a intermediate system X_1. The inertial system
%rotates arround z_I an Omega angle. 
%SECOND ROTAtion: The X_1 system is rotated (90º-I) arround x_1

L_1o=[1 0 0; 0 sin(I) cos(I); 0 -cos(I) sin(I)]; %1st rotation matrix. It transform a vector from local
%coordinates X_o to intemediate coordinates X_1. This matrix is the same
%for every orbital plane since all orbits have the same inclination

Re=6371;    %earth radius km. The distances are going to be refered to the earth radius
h=h/Re;     %high in earth radius
GM=6.474e-11*5.972e24*3600/Re^3; %Gravitational constant*earth mass [(Re^3)/min^2]
Ro=1+h;     %orbital radius in earth radius

w=sqrt(GM/(Ro^3));  %sat's angular speed [rad/min]
phi=0:2*pi/N_sat:(N_sat-1)*2*pi/N_sat; %relative angle between every sat of a plane respect the first one
f=8;                %factor to module the relative desfase of the first sat of one plane to the first sat of the next plane
l=length(t);

for p=1:N_planes    %going through each plane
    phi_0(p)=f*2*pi*(p-1)/(N_planes*N_sat);     %relative desfase of the first sat of one plane to the first sat of the next plane
    L_I1=[cos(Omega(p)) -sin(Omega(p)) 0; sin(Omega(p)) cos(Omega(p)) 0; 0 0 1]; %2nd rotation matrix. It transforms a vector from intermediate coordinates X_1
    %to inertial coordinates X_I. This matrix is characteristic of every
    %plane since evey plane has a diferent longitude of the ascending node
    %(Omega)
    L_eo=L_I1*L_1o;   %Global rotation matrix form local coordinates to nertial coordinates
    for s=1:N_sat   %going through every sat of the plane
        for k=1:l   %going through every time step
            X_loc=[Ro*cos(w*t(k)+phi(s)+phi_0(p));0;-Ro*sin(w*t(k)+phi(s)+phi_0(p))]; %Coordinates of the sat at the given time in local coordinates X_o.
            %It defines a circumference at the plane x_o-z_o with a angular
            %velocity w, a initial phase respect the first sat and a
            %initial phase respect the first sat of the previous orbit
            %first sat.
            X(:,(p-1)*N_sat+s,k)=L_eo*X_loc;    %Coordinates of the sat in inertial system (Coordinate [x y z], Sat number, instant(time)]
        end
    end        
end
end

