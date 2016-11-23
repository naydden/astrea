function tf = extra_lifetimeH(m0,S,cd,h0,hf)
    
% DESCRIPTION
% This function computes the time between altitudes h0 and hf
% hf is the design altitude. Includes postprocessing (Last lines)

% INPUT VARIABLES FORMAT
% m0 = Mass of the satellite assumed constant [kg]
% S  = Surface of the satellite [m2]
% cd = Drag coefficient.
% h0 = Starting height of the orbit [m]
% hf = Final height of the orbit [m]

% Test inputs
% Aplication example (Copy in command window)
% tf=extra_lifetimeH(4,0.1,2.2,600e3,500e3)

%% 1. Constants %%
%___Simulation Parameters____
h=1;   %Time step in terms of period (dt=h*P)

%___Physical Parameters____
u=3.986e14; %GM Earth
RE=6.378e6; %Earth Radius [m]

%% 2. initial Map
Rp=RE+h0;
t(1)=0; 
r(1)=Rp; %Satellite Orbital radius [m]
P(1)=2*pi*(r(1)^3/u)^.5; % Orbit Period [s]
dt=h*P(1);
 
%% 3. Solve the time steps
%Iterate satellite with time until it achieves the ground
i=1;
while r(i)>RE+hf
     i=i+1;
     
     rho=rho_formules(r(i-1));  
     dP=-3*pi*rho*r(i-1)*(cd*S/m0)*dt;  % Reduction in the period/rev
     
     %Compute the variation on the orbital period:
     P(i)=P(i-1)+dP;
     r(i)=(P(i)^2/(4*pi^2)*u)^(1/3); % Orbit Period [s]
     t(i)=t(i-1)+dt;
     dt=h*P(i);
     
     %show one timestep each 4000 iterations
%     if mod(i,4000)==0
%         fprintf('\nDías t=%g  h=%g km\n', t(i)/(86400),(r(i)-RE)/1000);
%     end
end

%% Function output in days
tf=t(i)/(86400); 

%% OPTIONAL - POSTPROCESSING
% fprintf('\nGround achieved after %g days\n', t(i)/86400);
% plot(t/86400,(r-RE)/1e3);
% grid;
% xlabel('time [days]');ylabel('altitude [km]');
% tit=['Satellite Orbital decay (m=' num2str(m0) ' kg, S_a_v_g=' num2str(S) ' m^2)']; 
% title(tit);

end