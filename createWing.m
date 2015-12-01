%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wing = createWing(top,bottom,SegNum,cordNum,filename,span,start,Controls,plane)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createWing: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the lattice structure for the calculation of
% the vortex lattice method. It convernts the finite element states
% into a usable lattice.
% NOTE: Current calculations are not for wings with sweep,
%       taper, or other variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	createMeanCord2015_8_10
%           MATLAB 5.2 std fcns
% Inputs:   plane - identifies which plane the wing lies on:
%               x-y=1, x-z=2,y-z=3

wing.span = span;
wing.SegNum = SegNum;
wing.cordNum = cordNum;
wing.L = span/SegNum*ones(wing.SegNum,1);
wing.meanCord = createMeanCord(filename,wing.cordNum);
wing.Controls = Controls;
wing.start = start;
wing.Theta = zeros(wing.SegNum+1,1);
wing.plane = plane;

switch plane
    case 1
        wing.Y = [0;cumsum(span/SegNum*ones(wing.SegNum,1))]-span/2;
        wing.cord = interp1(top(:,2),top(:,1),wing.Y,'pchip')-interp1(bottom(:,2),bottom(:,1),wing.Y,'pchip');
        wing.Z = zeros(wing.SegNum+1,1);
        wing.PhiZ = zeros(wing.SegNum+1,1);
        wing.X = interp1(bottom(:,2),bottom(:,1),wing.Y,'pchip');
        wing.PhiX = zeros(wing.SegNum+1,1);
    case 2
        wing.Z = [0;cumsum(wing.L)]+wing.start(3);
        wing.cord = abs(interp1(top(:,2),top(:,1),wing.Z,'pchip')-interp1(bottom(:,2),bottom(:,1),wing.Z,'pchip'));
        wing.Y = zeros(wing.SegNum+1,1);
        wing.PhiZ = zeros(wing.SegNum+1,1);
        wing.X = interp1(top(:,2),top(:,1),wing.Z,'pchip');
        wing.PhiX = zeros(wing.SegNum+1,1);
    case 3
        wing.X = [0;cumsum(wing.L)]+wing.start(1);
        wing.cord = interp1(top(:,2),top(:,1),wing.X,'pchip')-interp1(bottom(:,2),bottom(:,1),wing.X,'pchip');
        wing.PhiX = zeros(wing.SegNum+1,1);
        wing.Z = zeros(wing.SegNum+1,1);
        wing.PhiZ = zeros(wing.SegNum+1,1);
        wing.Y = zeros(wing.SegNum+1,1);
end