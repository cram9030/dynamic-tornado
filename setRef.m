%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ref]=setRef(L,bodyW,cord,bodyC,centroid,noseH)
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
% Inputs:   L - Length vector for wings
%           bodyW - Body Width
%           cord - Wing Cord Length
%           bodyC - Body Length
%           centroid - Center of mass
%           noseH - Length that the nose extends beyond the wing
% Output:   ref - Refrence structure containing:
%                   b_ref - Span
%                   S_ref - Refrence Area
%                   C_mgc - Mean Geometric Chord
%                   C_mac - Mean aerodymaic chord
%                   mac_pos - Mean aerodymaic chord position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the span
ref.b_ref = config('b_ref');
ref.b_ref = sum(L)+bodyW;

%Calculate the refrence area
ref.S_ref = config('S_ref');
ref.S_ref = sum(L)*cord+bodyW*bodyC;

%Calculate the Mean Geometric Cord
ref.C_mgc = ref.S_ref/ref.b_ref;

%Calculate the Mean Aerodynamic Cord
% EQN: \frac{2}{S} \int_0^{\frac{b}{2}c(y)^2dy
ref.C_mac = 2*(cord^2*sum(L)/2+bodyC^2*bodyW)/ref.S_ref;

%Calculate mac postion
% It is always be the centroid (or close to it as long as the body is thin
% compared to the wings) because there is no variance in the sections 
% cords Will need to be updated to accomidate future changes
ref.mac_pos = [-centroid(1)-noseH*2*(bodyC^2*bodyW)/ref.S_ref/ref.C_mac,centroid(2),0];
end