function [L] = Links(Xs,Xg,e_min)
%This function calculates the how many sats can a Ground Station see and 
%contact (link) in every instant

%Xs: sats coordinates in Earth radius in a inertial system centred at the center of the Earth. 
    %Xs(coordinates [x y z], sat number,instant(time))
%Xg: Ground Station coordinates in earth radius in a inertial system centred at the center of the Earth. . 
    %Xg(coordinates [x y z], instant(time))
%e_min: minimum elevation respect the horitzon from which the station can
%sees a sat. This elevaton angle resticts the view of a station because of
%the atmosphere and the obsacles. Measured in degrees

e_min=e_min*pi/180; %elevation angle in radians
l=size(Xg,2);       %number of instants to calculate
N_sat=size(Xs,2);   %number of satellites in the constellation

%for calulate the number of links ina given instant it will be computed the
%angle (alpha) between vector of position of the station and the vector
%station-sat.
alpha_lim=pi/2-e_min;  %maximum angle alpha for seeing a sat

L=zeros(1,l); %L will contain the number of links in every instant. It is
%defined as a array of zeros to fill in later

for i=1:l   %going through each instant 
    for j=1:N_sat  %for every sat
        Xgs=Xs(:,j,i)-Xg(:,i); %vector station-sat(j) at the instant i
        alpha=acos((Xgs'*Xg(:,i))/norm(Xgs)); %angle between the 
        %the vector position of the station and the vector station-sat(j)
        %at the instant i
        if alpha<alpha_lim %if the ange alpha is smaller of alpha_lim 
            L(i)=L(i)+1;    %the station can contact with the sat and it is
            %added a link at the instant i
        end
    end
end
end


