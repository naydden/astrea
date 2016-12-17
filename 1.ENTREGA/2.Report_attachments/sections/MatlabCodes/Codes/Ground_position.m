function [X] = Ground_position(lambda,long,t)
%Gorund_position calculates the positon of a place of the earth, with a
%wiven longitude and latitude, in an inertial system of coordiantes X_I in
%a given period of time

%lamda: is the latitude of the place in degrees[-90º,90º]
%long: is the longitude of the place iin degrees [0º,360º)
%t: is a array of time in minutes 

long=long*pi/180;   %longin radians
lambda=lambda*pi/180;   %lambda in radians

Rg=cos(lambda); %distance between the Earth point to the rotation axis of the Earth [0,1]
hg=sin(lambda); %distance between th Earth point to the equatorial plane [0,1]
%this 2 distances are measured in earth radius. 

w=2*pi/(24*60); %angular velocity of the Earth [rad/min]  

X=[Rg*cos(w*t+long);Rg*sin(w*t+long);hg*ones(1,length(t))]; %position of the point in X_I
%for the given time. The point descrives a cicumference with a Rg radius in the horizontal
%plane x_I-y_I with an angular velocity w and a initial phase long. The
%coordinate z_I is allways hg. It is expresed in earth radius. 
%X(cooridinates, instant (time)

end

