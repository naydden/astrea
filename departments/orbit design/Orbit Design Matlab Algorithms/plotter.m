%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Routine created to plot the results for several altitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
close all;
clc;
%%
hatm = [113 200 300];  %[km]

N = 20;
phi=zeros(length(hatm),N);
d=phi;
h=zeros(1,N);

for j=1:length(hatm)
    for i=1:N
        h(i) = 400+25*i;
        [phi(j,i), d(j,i)] = satsatVisibility(hatm(j),h(i));
    end
end

%% Post-processing
figure;
ylabel('angle between sattellites [^{\circ}]'); xlabel('altitude [km]');
grid;
hold on;
s=[];
for i=1:length(hatm)
    plot(h,radtodeg(phi(i,:)));
    s=[s; 'h_a_t_m = ' num2str(hatm(i))];
end
legend(s);