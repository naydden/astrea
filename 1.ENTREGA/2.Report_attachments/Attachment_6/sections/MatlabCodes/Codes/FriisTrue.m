%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelors Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% SATELLITE DESIGN TEAM
% FRIIS EQUATION RANGE CALCULATOR

%{
This matlab routine calculates the strength of the received signal for
distance variation. Taking into account losses considered on the Link
Budget Annex. Calculus are based on Friis equation and Bandwidth-Strength
correlation.
%}
clc; clear; close all;
%
%% Initial variables definition
r_start = 100; %km
r_end = 1000; %km
n=1; %step
Lprop=0;

%% Sband Space to Ground station calculus
Gt=10; %dBi
Gr=35.4; %dBi
f=2228; %MHz
Pw=2;
L_cables = 2; Labs=0.58; Laml=1; Lph=20;
Lprop = L_cables+Labs+Laml+Lph;
[ PrS,radiusS, FSL_S] = friis_eq( Gt, Gr, f, Pw,r_start, r_end, n,Lprop);

%% Sband intersatellite
Gt=10; %dBi
Gr=10; %dBi
f=2228; %MHz
Pw=2;
%L_cables = 2; Labs=0.58; Laml=1; Lph=20;
Lprop = 0;%L_cables+Labs+Laml+Lph;
[ PrSi,radiusSi, FSL_Si] = friis_eq( Gt, Gr, f, Pw,r_start, r_end, n,Lprop);

%% Xband signal strength between satellites calculation
Gt=16.5; %dB
Gr=16.5; %dB
f=9e3; %MHz
Pw=12; %W
[ PrX1,radiusX1, FSL_X1] = friis_eq( Gt, Gr, f, Pw,r_start, r_end, n,Lprop);

%% Xband between sat to ground
r_start = 100; %km
r_end = 1000; %km
n=1; %step
Gt=16.5; %dB
Gr=36; %dB
f=10e3; %MHz
Pw=12; %W
L_cables = 2; Labs=0.58; Laml=1; Lph=20;
Lprop = L_cables+Labs+Laml+Lph;
[ PrX2,radiusX2, FSL_X2] = friis_eq( Gt, Gr, f, Pw,r_start, r_end, n, Lprop);

%% PLOT results
X1=radiusS;
YMatrix1=[PrS' PrSi' PrX1' PrX2'];
llegenda1='Sband Sat to GND';
llegenda2='Sband Sat to Sat';
llegenda3='Xband Sat to Sat ';
llegenda4='Xband Sat to GND';
eixX='Distance [Km]';
eixY='Power [dB]';
titol='Change in received power';
createfigure(X1, YMatrix1,llegenda1,llegenda2,llegenda3,llegenda4, titol,eixX, eixY)
legend(llegenda1,llegenda2,llegenda3,llegenda4)
