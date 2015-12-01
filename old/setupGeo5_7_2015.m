%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [geo] = setupGeo5_7_2015(ref,CG,span,cordNum,SegNum,numSpanB)

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
% Calls:	
% Inputs:   ref - refrence point
%           CG - center of gravity
% Output:   geo - geomotry structure containing:
%                   ref - refrence point
%                   CG - center of gravity
%                   b - total span of the wing
%                   nx - number of cord wise panels
%                   ny - number of spanwise panels
%                   nwing - number of wings, set to 1 will need to be
%                           updated in the future
%                   startx - begining x coordinate
%                   starty - begining y coordinate
%                   startz - begining x coordinate
%                   semetric - wing is always assumed to not be symmetric
%                   numSpanW - number of wing span panels
%                   numSpanB - number of body span panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geo.ref_point = ref;
geo.CG = CG;
geo.b = span;
geo.nx = cordNum;
geo.ny = SegNum+numSpanB;
geo.numSpanW = SegNum;
geo.numSpanB = numSpanB;
geo.nwing = 1;
geo.startx = 0;
geo.starty = 0;
geo.startz = 0;
geo.symetric = 0;