%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ref]=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setRef: Essential function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the refrence values of the wing				
% NOTE: Current calculations are not for wings with sweep,
%       taper, or other variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	generateLattice5_6												
% Calls:	MATLAB 5.2 std fcns
% Inputs:   S_ref - refrence area
%           C_ref - refrence chord
%           B_ref - refrence span
% Output:   ref - Refrence structure containing:
%                   b_ref - Span
%                   S_ref - Refrence Area
%                   C_mgc - Mean Geometric Chord
%                   C_mac - Mean aerodymaic chord
%                   mac_pos - Mean aerodymaic chord position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the span
ref.b_ref = config('b_ref');
ref.b_ref = B_ref;

%Calculate the refrence area
ref.S_ref = config('S_ref');
ref.S_ref = S_ref;

%Calculate the Mean Geometric Cord
ref.C_mgc = ref.S_ref/ref.b_ref;

%Calculate the Mean Aerodynamic Cord
ref.C_mac = C_ref;

%Calculate mac postion
% It is always be the centroid (or close to it as long as the body is thin
% compared to the wings) because there is no variance in the sections 
% cords Will need to be updated to accomidate future changes
ref.mac_pos = mac_pos;
end