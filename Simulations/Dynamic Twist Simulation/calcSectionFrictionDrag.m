%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fw = calcSectionFrictionDrag(U,t,c,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcFrictionDrag: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the skin friction coefficent of drag of the wind
% tunnel model of MADCAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
% Inputs:   U - wind speed
%           t - thickness
%           c - cord
%           s - spanwise section
% Output:   Cd - Coefficient of drag from skin friction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Equations from Fundamentals of Airplane Flight Mechanics and Aerodynamics
%for Engineers

%This configuration does not seem to work. The combination of u,s,and c do
%not seem to allow it to be simply added together and it is not worth it to
%figure out the appropriate translation right now

%Body and wing deminsions
Sexposed = 2*s*c;

%Refrence Area
S = 1.289685*0.3048;

%Air constants
airDensity = 1.225; %kg/m^3
airDynamicViscosity = 12.27e-6; %Pa*s
airSoundSpeed = 340.29; %m/s
mach = U/airSoundSpeed;

%Reynolds number for body and wings
ReWings = airDensity*U*c/airDynamicViscosity;

%Compressibility factor
CF_k = (1.0+0.2*mach.^2).^(-0.467);

%Interferance Factor
Iw = 1.20;

%Form Factor
FFw = 1.0+1.6*(t/c)+100*(t/c)^4;

%Wetted area
Swing = 2.0*(1+0.2*t/c)*Sexposed;

%Coefficient of Friction - Prandtl-Schlichting formula transition region
Cfw = 0.455./(log10(ReWings)).^2.58-1700./ReWings;

%Coefficient of Friction - Prandtl-Schlichting formula
Cfw = 0.455./(log10(ReWings)).^2.58;

%Force Calculation
fw = Cfw.*CF_k*Iw*FFw*Swing;

%CD Calculation
Cd = 1.1*fw/S;