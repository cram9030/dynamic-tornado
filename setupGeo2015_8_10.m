%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [geo] = setupGeo2015_8_10(ref_point,CG,Wings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setupState5_7_2015: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the geo structure for the geomotry
% NOTE: The traditional version of TORNADO has all of the wing geomotry
% contained with in the geo structure. In Dynamic TORNADO that information
% is soley contained in the Lattice and the Geo is only used for refrence
% points and CG but might need to be updated latter to work correctly with
% other parts of pre-existing code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	createMeanCord2015_8_10
% Inputs:   ref - refrence point
%           CG - center of gravity
%           nwing - number of wings
%           Wings - Array of wing structures
%                   wing - Wing structure
%                           meanCord - mean cord line
%                           cord - cord length
%                           span - wing span
%                           L - Length array
%                           SegNum - number of spanwise segments
%                           cordNum - number of cordwise segments
%                           Theta - wing twist array
%                           Z - z displacement array
%                           PhiZ - Z slope array
%                           X - x displacement array
%                           PhiX - x slope array
%                           start - wing start location
% Output:   geo - Geometry of the plane
%                   ref_point - aircraft geometry refrence point
%                   CG - aircraft center of gravity
%                   nwing - number of wings
%                   symetric - is wing symetric about y axis
%                   Wings - Array of wing structures
%                           wing - Wing structure
%                                   meanCord - mean cord line
%                                   cord - cord length
%                                   span - wing span
%                                   L - Length array
%                                   SegNum - number of spanwise segments
%                                   cordNum - number of cordwise segments
%                                   Theta - wing twist array
%                                   Z - z displacement array
%                                   PhiZ - Z slope array
%                                   X - x displacement array
%                                   PhiX - x slope array
%                                   start - wing start location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geo.ref_point = ref_point;
geo.CG = CG;
geo.nwing = length(Wings);
geo.symetric = 0;
geo.Wings = Wings;
geo.nx = 0;
geo.ny = 0;
for i = 1:length(Wings)
    geo.nx = Wings(i).wing.cordNum+geo.nx;
    geo.ny = Wings(i).wing.SegNum+geo.ny;
end