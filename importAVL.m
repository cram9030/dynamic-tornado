%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state,geo,lattice,ref] = importAVL(fileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importAVL: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads and AVL file and generates dynamic tornado state,
% geomotry, lattice, and refrence strcutures that can be used for
% simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	
% Inputs:   fileName - path and name of AVL file
%           
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    fid = fopen(fname);
catch exception
    throw(exception)
end

currentLine = tline(fid);
lineNum = 1;

while currentLine ~= -1
    switch lineNum
        case 1
        case 2
            U_inf = str2num(currentLine)*340.29;
        case 3
            lineSplit = strsplit(currentLine);
            ref_point(1) = str2num(lineSplit{1});
            ref_point(2) = str2num(lineSplit{2});
            ref_point(3) = str2num(lineSplit{3});
        case 4
            lineSplit = strsplit(currentLine);
            S_ref = str2num(lineSplit{1});
            C_ref = str2num(lineSplit{2});
            B_ref = str2num(lineSplit{3});
        case 5
            lineSplit = strsplit(currentLine);
            CG(1) = str2num(lineSplit{1});
            CG(2) = str2num(lineSplit{2});
            CG(3) = str2num(lineSplit{3});
        case 6
            Cdp = str2num(currentLine);
        otherwise
            if strcmp(currentLine,'SURFACE')
                surfaceLine = tline(fid);
                surfaceLineNum = 1;
                while ~strcmp(surfaceLine, 'SURFACE');
                    switch surfaceLineNum
                        case 1
                        case 2
                        case 3
                            lineSplit = strsplit(currentLine);
                            SegNum = str2num(lineSplit{1});
                            
                    surfaceLine = tline(fid);
    end
    lineNum = lineNum + 1;
    currentLine = tline(fid);
end

fclose(fid)