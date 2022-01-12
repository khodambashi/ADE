%% Introduction
% This code solve the classical Advection-Dispersion Equation (ADE) for uniform
% channel and flow  parameters with pdepe function in MATLAB.

% This code available at :
%                      https://github.com/khodambashi/ADE.git

% For more information : 
%                      https://www.mathworks.com/help/matlab/ref/pdepe.html

% sajad Khodambashi Emami, Ver 1.0, jan 2022

%% Clear Command Window and Workespace, Close figuers 
close all
clear
clc

%% Difine Models Parameters
global D U S K

D = 0.24; % Velocity [m/s]
U = 0.37; % Longitudinal Dispersion Cofficent [m^2/s]
S = 0;    %  Source term
K = 0 ;   %  First order decay coefficient 

%% Spatial Mesh

Xmin = 0 ;                % pollutant injection site
L = 281 ;                 % Length of river [m]
NX = 282 ;                % Number of spatial subdivisions
X = linspace(Xmin,L,NX) ; % Subdivision of the Spatial domain.
dX = L/(NX-1) ;           % space-step

%% Temporal Mesh

Tmin = 0 ;                   % initial time
Tmax = (L/U)*1.5  ;          % lengh of simulation [s]
NT = 201 ;                   % Number of temporal subdivisions   
T = linspace(Tmin,Tmax,NT) ; % Subdivision of the temporal domain.
dT = Tmax/(NT-1) ;           % time-step


%% Time Pattern of Pollutant Injection
global t c
c = [10 ,10 ,0 ,0];
t = [0 ,100 ,150 ,10000];

%%  Symmetry constant
m = 0 ;

%% Main Solve Function
C = pdepe(m,@ADEfunc,@ADEic,@ADEbc,X,T);

%% plot plluutant transport in river

for i=1:NT
    
   plot(X,C(i,:),'k','LineWidth',1) 
   title('Pllutant Transport in River')
   xlabel('River Route [m]')
   ylabel ('Concentration [ppm]')
   box on
   grid on
   
    axis([0 L 0 max(c)])
    pause(0.0412)
end
close all
%% Difine Equation for pdepe
function [c,f,s] = ADEfunc(~,~,u,dudx)
global D U S K

   c = 1;
    
   f = D * dudx;
    
   s = -U * dudx - K * u + S;

end

%% Initial Condition
function u0 = ADEic(~)

   u0 = 0;

end

%% Boundray Condition
function [pl,ql,pr,qr] = ADEbc(~,ul,~,~,T)
global t c D
C_in = interp1(t,c,T) ;

   % Left Boundary Conditions (time pattern of pollutant injection)
 pl = ul - C_in  ; 
 ql = 0 ; 
 
   % Right Boundary Conditions ()
 pr = 0 ;
 qr = 1/D ; 
 
end

%%-------------------------------------------------------------------------