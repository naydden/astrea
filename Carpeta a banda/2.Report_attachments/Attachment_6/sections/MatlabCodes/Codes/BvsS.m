%ASTREA CONSTELLATION
%PROJECTS - 220028
%Aerospace Engineering Barchelors Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% SATELLITE DESIGN TEAM
% SENSITIVITY vs. BANDWITH

%{
This matlab routine calculates the Bandwidth to received signal Strength
correlation.
%}
clc; clear; close all;
% 
%% Sband calculus
C = 25*8*1e6 ; %bits/s
B = 100; %Hz
N = 1.38e-23*87.79*B; %W
Binf=5e6; Bsup=56e6; Bstep=1e4;

i=1;
S=zeros(size(Binf:Bstep:Bsup));
for B = Binf:Bstep:Bsup
    N = 1.38e-23*87.79*B;
    S(i) = 10*log10(100*(N*(2^(C/B)-1)));
    %SdB(i) = 10*log10*(1000*S/1);
    i=i+1;
end
Bs = Binf:Bstep:Bsup;
Ss = S;

%% Xband calculus
Binf=5e6; Bsup=100e6; Bstep=1e4;

i=1;
S=zeros(size(Binf:Bstep:Bsup));
for B = Binf:Bstep:Bsup
    N = 1.38e-23*87.79*B;
    S(i) = 10*log10(100*(N*(2^(C/B)-1)));
    %SdB(i) = 10*log10*(1000*S/1);
    i=i+1;
end
Sx=S; Bx = Binf:Bstep:Bsup;
%% PLOT results
singleFig(Bx, Sx);
