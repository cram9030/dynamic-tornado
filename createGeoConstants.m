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
% Note: Accending/decending ordering is important when using trapz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CA = abs(trapz(AF.UpperX,AF.UpperY)+trapz(AF.LowerX,AF.LowerY));
centroid = 1/CA*[-trapz(AF.UpperX,AF.UpperX.*AF.UpperY)-trapz(AF.LowerX,AF.LowerX.*AF.LowerY),-trapz(AF.UpperX,AF.UpperY.^2)-trapz(AF.LowerX,AF.LowerY.^2)];
I_x = (-trapz(AF.UpperX,AF.UpperY.*(AF.UpperX-centroid(1)).^2)-trapz(AF.LowerX,AF.LowerY.*(AF.LowerX-centroid(1)).^2))+CA*(ActLoc(1)-centroid(1))^2;
LowerY = interp1(AF.LowerX,AF.LowerY,AF.UpperX,'pchip');
I_z = trapz(AF.UpperY,AF.UpperX.*AF.UpperY.^2)-2*trapz(AF.UpperY,AF.UpperX.*LowerY.*AF.UpperX)+trapz(AF.LowerY,AF.LowerX.*AF.LowerY.^2)+CA*(ActLoc(2)-centroid(2))^2;
J_z = I_x+I_z;