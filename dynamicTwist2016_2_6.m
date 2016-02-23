%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,G] = dynamicTwist2016_2_6(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamicTwist2016_2_6: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates the dynamic tip twist of the wing for a
% proscribed tip profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
%           seitupState5_7_2015
%           dynamicSolver6_23_2015
%           setupGeo6_15_2013
%           setRef6_15
% Inputs:  
% Output:   CL - array lift coefficent at the times from twist input
%           CD - array drag coefficent at the times from twist input
%           LD - array lift/drag ratio at the times from twist input
%           forceTwist - array of required torque for specified tip twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tunnelWing = [];
airSpeed = [];
cordNum = [];
ActLoc = [];
C_ref = [];
rhom_inf = [];
centroid = [];
S_ref = [];
B_ref = [];
alpha_root = [];
sideSlip = [];

load('OptimizationLoadFile2016_2_21.mat')

[CL,CD,Cm,Cl,Cn,LD] = dynamicTwistOutput2016_2_23(x);

dt = 1/50;
twist(:,1) = x(1:length(x));
twistRate(:,1) = diff(twist(:,1))/dt;

F = [mean(CD*.5*rhom_inf*U_inf^2*S_ref);
    mean(CL)*.5*rhom_inf*U_inf^2*S_ref;
    max(abs(twist(:,1)))/dt;
    x(1)-x(end);
    CL*.5*rhom_inf*U_inf^2*S_ref;
    ];

mean(CD*.5*rhom_inf*U_inf^2*S_ref)

% Define the derivatives.
G=[];