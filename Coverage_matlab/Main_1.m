%%
clc; clear;
N_sat=21;
N_planes=14;
h=545;
I=18;
phase=360/(N_planes-1);
t=0:5:24*60*2;
e_min=20;
lambda=0:30:90; lat=length(lambda);
mu=0; long=length(mu);
l=length(t);

Xs=Orbital_position(N_sat,N_planes,h,I,phase,t);
L=zeros(lat,l);

for i=1:lat
    Xg=Ground_position(lambda(i),mu,t);
    [L(i,:)]=Links(Xs,Xg,e_min);
end

%%
figure; plot(t/60,L); legend 0 30 60 90 ; xlabel('time (h)'); ylabel('links');

