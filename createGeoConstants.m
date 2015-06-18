%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [centroid,I_z,I_x,J_z,CA] = createGeoConstants(AF,ActLoc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generateLattice5_6: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the centroid, area moments of interia,
% and the cross sectional area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
% Inputs:   AF - Airfoil structure containing:
%               UpperX - X position for top of airfoil
%               UpperY - Y position for top of airfoil
%               LowerX - X position for the bottom of airfoil
%               LowerY - Y position for the bottom of airfoil
%               Name - Name of airfoil
%           ActLoc - Actuator location
% Output:   centroid - centroid of airfoil
%           I_z - z area moment of inertia 
%           I_x - x area moment of inertia 
%           J_z - combined rotational z area moment of inertia 
%           CA - cross sectional area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CA = trapz(AF.UpperX,AF.UpperY-AF.LowerY);
centroid = [1/CA*trapz(AF.UpperX,AF.UpperX.*(AF.UpperY-AF.LowerY)),0];
I_x = trapz(AF.UpperX,(AF.UpperY-AF.LowerY).*(AF.UpperX-centroid(1)).^2)+CA*(ActLoc(1)-centroid(1))^2;
I_z = 2*trapz(AF.LowerY,AF.UpperX.*AF.UpperY.^2)+CA*(ActLoc(2)-centroid(2))^2;
J_z = I_x+I_z;