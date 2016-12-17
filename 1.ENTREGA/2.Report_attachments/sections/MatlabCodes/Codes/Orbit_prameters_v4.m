%----ASTREA CONSTELLATION----
%PROJECTS - 220028
%Aerospace Engineering Barchelor's Degree
%ESEIAAT - UPC
%Autumn 2016-2017

% ORBIT DESIGN TEAM
% ORBIT PARAMETERS

% This routine calculates the minimum number of plans and satellites
% to guarantee global coverage.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% _v2 --> Computes angles with arccosine
% _v3 --> Also iterates in f

clc
clear all
close all

tic

%% DATA INPUT

Re = 6371;              %Earth's Radius [Km]
Se = 4*pi*Re^2;         %Earth's Surface [Km^2]
h = 500:5:600;          %Height vector [km]
in = 75:3:90;           %Inclination vector [�]
pmin = 5; p=12;         %Minimum number of plans
pmax = 20;              %Maximum number of plans
sppmin = 10; spp=10;    %Minimum number of satellites each plane
sppmax = 24;            %Maximum number of satellites each plane
typeWD = 2;             %Defines:
                        % 1- Semi-Walker-Delta
                        % 2- 3/4 of Walker-Delta
                        % 3- Full Walker-Delta
plaunch = 5.76e6;       %Launch price
psat = 250e3;           %Price per satellite
eo = 20;                %Vison angle
q=0;                    %Variable used as counter
nplanes_win=0;          %Boolean, Indicates a succesful combination with that number of planes

%% COMPUTATIONS
%result=zeros(6,100);

for j=1:length(h)
    fprintf('\n\nH = %g km',h(j))
    Bo = visibleangle(h(j),eo);  %Footprint angle
                    
    for k=1:length(in)
        fprintf('\ni=%g�:  ',in(k))
        
        for p=pmin:pmax
            F=0:1:p-1; nF=length(F); anglesp_f=zeros(1,nF);
            
            for spp=sppmin:sppmax
                angles = 360/spp; %Angle between satellites of same plane
                
              for kf=1:nF
               % We look for the optimum f for this spp,npp,inc,h
               anglesp_f(kf)=distance_opt_f(in(k),h(j),p,spp,typeWD,F(kf));   
               %Angle between satellites of different planes
              end
                
                [anglesp,kf_min]=min(anglesp_f);

                 if (angles <= 2*Bo) && (anglesp <= Bo)
                    
                    q=q+1;
                    result(:,q)=[h(j) in(k) spp p spp*p F(kf_min)];
                    fprintf(['Solution %g: h=%gkm - i=%g� - '...
                       '%gspp - %gplanes - %g sats - f=%g'],q,result(:,q))
                    nplanes_win=1;
                    break %Stop iteration in satellites per plain
       
                 end
                 
            end
            if nplanes_win     %If a combination with that num of planes
                nplanes_win=0; %works, then jump to the next inclination
                break
            end
        end
    end
end

%% Optimum Results
% Minimum number of cubesats
numsat = result (5,1:q);
[minsat,index] = min(numsat);

disp(sprintf('\nMINIMUM SATELLITES\nheight: %d', result(1,index)))
disp(sprintf('inclination: %d',result(2,index)))
disp(sprintf('number of satellites in each plane: %d',result(3,index)))
disp(sprintf('number of planes: %d',result(4,index)))
disp(sprintf('minimum number of satellites: %d',result(5,index)))

% Minimum number of planes
numplane = result(4,1:q);
[minplane,index]= min(numplane);

disp(sprintf('\nMINIMUM PLANES\nheight: %d', result(1,index)))
disp(sprintf('inclination: %d',result(2,index)))
disp(sprintf('number of satellites in each plane: %d',result(3,index)))
disp(sprintf('number of planes: %d',result(4,index)))
disp(sprintf('minimum number of satellites: %d',result(5,index)))

% Price Optimization
TCost = plaunch*result(4,:) + psat*result(5,:); TCost = TCost/(1e6);
[CostOpt,index] = min(TCost);
height = result(1,:);

disp(sprintf('\nMINIMUM PRICE\nprice[M]: %d', TCost(index)))
disp(sprintf('height: %d', result(1,index)))
disp(sprintf('inclination: %d',result(2,index)))
disp(sprintf('number of satellites in each plane: %d',result(3,index)))
disp(sprintf('number of planes: %d',result(4,index)))
disp(sprintf('minimum number of satellites: %d',result(5,index)))
        
fprintf('\n\nCAREFUL!! MORE THAN ONE SOLUTION COULD HAVE THE SAME MIN NUMBER\n')     

%% Plots

x=0:q-1;
plot(x,numplane,x,numsat,x,TCost,x,height)
strValues = strtrim(cellstr(num2str([x(:) numplane(:)],'(%d, %d)')));
text(x,numplane,strValues,'VerticalAlignment','bottom');
strValues = strtrim(cellstr(num2str([x(:) numsat(:)],'(%d, %d)')));
text(x,numsat,strValues,'VerticalAlignment','bottom');
strValues = strtrim(cellstr(num2str([x(:) TCost(:)],'(%d, %.1f)')));
text(x,TCost,strValues,'VerticalAlignment','bottom');
strValues = strtrim(cellstr(num2str([x(:) height(:)],'(%d, %d)')));
text(x,height,strValues,'VerticalAlignment','bottom');
legend('number of planes','number of satellites','Total Cost (M)')



toc