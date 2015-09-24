%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NP] = findVerticalNeutralPoint(delta,state,geo,lattice,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findNeutralPoint: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the neutral point of the aircraft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	findMDerv
% Inputs:   delta - alpha difference for simulations
%           geo - Geometry of the plane
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
%           state - state structure containing:
%                   alpha - local angle of attack of wings
%                   alpha_root - angle of attack at wing root and boyd
%                   betha - angle of side sweep
%                   P - Roll angular velocity [deg/s], set to zero
%                   Q - Pitch angular velocity [deg/s], set to zero
%                   R - Yaw angular velocity [deg/s], set to zero
%                   adot - Angle of attack time derivative, (Alpha_dot),
%                   [deg/s], set to zero
%                   bdot - Angle of sideslip time derivative, (Beta_dot),
%                   [deg/s], set to zero
%                   AS - local Airspeed
%                   U_inf - fair field airspeed
%                   rho - Air density
%                   ALT - altitude, set to zero because air density is
%                   provided
%                   alphadot - likely legacy code, set to zero
%                   bethadot - likely legacy code, set to zero
%                   pgcorr - Apply Prandtl-Glauert Correction [0 1], zero
%                   is false, one is true, set to zero
%           lattice - lattice structure containing:
%                   COLLOC - cordNum*(SegNum+numSpanB)by 3 positional
%                   matrix of collocation points
%                   VORTEX - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                   containing x,y,z positions for each corner point of the
%                   vortex
%                   N - Normal vectors of VORTEX
%                   XYZ - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                   containing x,y,z positions for each external corner
%           ref - Refrence structure containing:
%                   b_ref - Span
%                   S_ref - Refrence Area
%                   C_mgc - Mean Geometric Chord
%                   C_mac - Mean aerodymaic chord
%                   mac_pos - Mean aerodymaic chord position
% Output:   neutral point - cordwise neutral point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a.delta = delta;
a.state = state;
a.geo = geo;
a.lattice = lattice;
a.ref = ref;

NP = fminsearch(@(x) findMDerv(x,a),[-0.01,0]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findMDerv: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the moment derivative for neutral point search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	aeroDerivative
% Inputs:   x - cordwise displacement
%           a - temp variable containing required content
%               delta - alpha difference for simulations
%               geo - Geometry of the plane
%                       ref_point - aircraft geometry refrence point
%                       CG - aircraft center of gravity
%                       nwing - number of wings
%                       symetric - is wing symetric about y axis
%                       Wings - Array of wing structures
%                               wing - Wing structure
%                                       meanCord - mean cord line
%                                       cord - cord length
%                                       span - wing span
%                                       L - Length array
%                                       SegNum - number of spanwise segments
%                                       cordNum - number of cordwise segments
%                                       Theta - wing twist array
%                                       Z - z displacement array
%                                       PhiZ - Z slope array
%                                       X - x displacement array
%                                       PhiX - x slope array
%                                       start - wing start location
%               state - state structure containing:
%                       alpha - local angle of attack of wings
%                       alpha_root - angle of attack at wing root and boyd
%                       betha - angle of side sweep
%                       P - Roll angular velocity [deg/s], set to zero
%                       Q - Pitch angular velocity [deg/s], set to zero
%                       R - Yaw angular velocity [deg/s], set to zero
%                       adot - Angle of attack time derivative, (Alpha_dot),
%                       [deg/s], set to zero
%                       bdot - Angle of sideslip time derivative, (Beta_dot),
%                       [deg/s], set to zero
%                       AS - local Airspeed
%                       U_inf - fair field airspeed
%                       rho - Air density
%                       ALT - altitude, set to zero because air density is
%                       provided
%                       alphadot - likely legacy code, set to zero
%                       bethadot - likely legacy code, set to zero
%                       pgcorr - Apply Prandtl-Glauert Correction [0 1], zero
%                       is false, one is true, set to zero
%               lattice - lattice structure containing:
%                       COLLOC - cordNum*(SegNum+numSpanB)by 3 positional
%                       matrix of collocation points
%                       VORTEX - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                       containing x,y,z positions for each corner point of the
%                       vortex
%                       N - Normal vectors of VORTEX
%                       XYZ - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                       containing x,y,z positions for each external corner
%               ref - Refrence structure containing:
%                       b_ref - Span
%                       S_ref - Refrence Area
%                       C_mgc - Mean Geometric Chord
%                       C_mac - Mean aerodymaic chord
%                       mac_pos - Mean aerodymaic chord position
% Output:   Cm_a - moment derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cm_a = findMDerv(x,a)

delta = a.delta;
state = a.state;
geo = a.geo;
lattice = a.lattice;
ref = a.ref;
state.alpha_root = x(2);

[geo] = setupGeo2015_8_10(geo.ref_point,[geo.CG(1),0,x(1)],geo.nwing,geo.Wings);
[results] = aeroDerivative(delta,state,geo,lattice,ref);
Cm_a = abs(results.Cm_a);
disp(['CG: ',num2str(x(1))])
disp(['AoA: ',num2str(x(2))])
disp(['Cm_a: ',num2str(results.Cm_a)])

end