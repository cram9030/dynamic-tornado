%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results] = aeroDerivative(delta,state,geo,lattice,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aeroDerivative: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the aerodynamic derivaties
% NOTE: Currently it is only set up to do derivatives with refrence to
% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	dynamicSolver6_23_2015
%           coeff_create2015_8_13
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
% Output:   results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_range = [state.alpha_root-delta,state.alpha_root+delta];
[results]= dynamicSolver6_23_2015(state,geo,lattice);

for i = 1:length(alpha_range)
    [state] = setupState2015_8_10(alpha_range(i)*ones(size(state.alpha)),alpha_range(i),state.betha,state.AS,state.U_inf,state.rho);
    [lattice] = generateLattice2015_8_10(geo,0,0,0,0,state);
    [tempResults]= dynamicSolver6_23_2015(state,geo,lattice);
    results.FORCE(i+1,:,:) = tempResults.FORCE;
    results.MOMENTS(i+1,:,:) = tempResults.MOMENTS;
end

[results]=coeff_create2015_8_17(results,lattice,state,ref,geo,delta);