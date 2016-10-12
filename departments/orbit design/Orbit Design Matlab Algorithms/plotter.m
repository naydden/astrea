%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine created to plot the results for several altitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
close all;
clc;
%%
hatm = 100;  %[km]

N = 20;
phi=zeros(1,N);
d=phi;
h=phi;

for i=1:N
    h(i) = 400+25*i;
    [phi(i), d(i)] = satsatVisibility(hatm,h(i));
end

%% Post-processing
plot(h,radtodeg(phi));
ylabel('angle between sattellites [^{\circ}]'); xlabel('altitude [km]');
grid;

    