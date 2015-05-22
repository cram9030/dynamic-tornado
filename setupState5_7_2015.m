%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state] = setupState5_7_2015(alpha,beta,airSpeed,airDensity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setupState5_7_2015: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the state structure for the environmental
% conditions for the simulation
% NOTE: All the dynamic states are set to zero and this function will need
% to be altered for them to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
% Inputs:   alpha - angle of attack, degrees
%           beta - angle of sidesweep, degrees
%           airSpeed - Airspeed in m/s
%           airDensity = Air density in kg/m^3
% Output:   state - state structure containing:
%                   alpha - angle of attack
%                   betha - angle of side sweep
%                   P - Roll angular velocity [deg/s], set to zero
%                   Q - Pitch angular velocity [deg/s], set to zero
%                   R - Yaw angular velocity [deg/s], set to zero
%                   adot - Angle of attack time derivative, (Alpha_dot),
%                   [deg/s], set to zero
%                   bdot - Angle of sideslip time derivative, (Beta_dot),
%                   [deg/s], set to zero
%                   AS - Airspeed
%                   rho - Air density
%                   ALT - altitude, set to zero because air density is
%                   provided
%                   alphadot - likely legacy code, set to zero
%                   bethadot - likely legacy code, set to zero
%                   pgcorr - Apply Prandtl-Glauert Correction [0 1], zero
%                   is false, one is true, set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state.alpha = alpha*pi/180;
state.betha = beta*pi/180;
state.P = 0;
state.Q = 0;
state.R = 0;
state.adot = 0;
state.bdot = 0;
state.AS = airSpeed;
state.rho = airDensity;
state.ALT = 0;
state.alphadot = 0;
state.bethadot = 0;
state.pgcorr = 0;