%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cd = calcFrictionDrag(U)
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
% Output:   Cd - Coefficient of drag from skin friction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Equations from Fundamentals of Airplane Flight Mechanics and Aerodynamics
%for Engineers

%Body and wing deminsions
t = 36.576e-3;
c = 0.3048;
Sexposed = 2*0.58205878*c;
crossSecWing = 0.007328002;
l = 0.4191;
d = 0.10374122;
h = 0.127;
noseSA = 0.0378896836;

%Refrence Area
S = 1.289685*c;

%Air constants
airDensity = 1.225; %kg/m^3
airDynamicViscosity = 12.27e-6; %Pa*s
airSoundSpeed = 340.29; %m/s
mach = U/airSoundSpeed;

%Reynolds number for body and wings
ReBody = airDensity*U*.4191/airDynamicViscosity;
ReWings = airDensity*U*0.3048/airDynamicViscosity;

%Compressibility factor
CF_k = (1.0+0.2*mach.^2).^(-0.467);

%Interferance Factor
Ib = 1.20;
Iw = 1.20;

%Form Factor
FFw = 1.0+1.6*(t/c)+100*(t/c)^4;
FFB = 1.0+60/(l/d)^3+0.0025*(l/d);

%Wetted area
Swing = 2.0*(1+0.2*t/c)*Sexposed;
SBody = noseSA-2*crossSecWing+2*l*d+2*h*l;

%Coefficient of Friction - Prandtl-Schlichting formula
Cfw = 0.455./(log10(ReWings)).^2.58-1700./ReWings;
Cfb = 0.455./(log10(ReBody)).^2.58-1700./ReBody;

%Force Calculation
fw = Cfw.*CF_k*Iw*FFw*Swing;
fb = Cfb.*CF_k*Ib*FFB*SBody;

%CD Calculation
Cd = 1.1*(fw+fb)/S;