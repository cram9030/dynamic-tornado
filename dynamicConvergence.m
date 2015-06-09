%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fa,States,DT] = dynamicConvergence(Ain,Bin,Fcin,SegNum,States,convergePercent,cord,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,state,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamicConvergence: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gives the aerodynamic forces and beam states
% that result from convergence of the vortex lattice method
% and the galerkin finite element method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 14a std fcns
%           generateLattice5_30
%           setupGeo5_7_2015
%           solver5_7_2015
%           coeff_create5_7_2015
% Inputs:   A - Galkerkin beam states matrix
%           B - Galkerkin beam input matrix
%           Fc - Control force
%           States - Current wing states
%           convergePercent - required percentage convergence
%           cord - wing cord
%           centroid - wing centroid
%           cordNum - number of segments that cords are split into
%           L - array of segments lengths
%           bw - body width
%           bc - body cord
%           noseH - nose cone distance
%           zb - body offset
%           numSpanB - number of body span sections
%           state - aerodynamic state
%           dt - maximum time step
% Output:   Fa - Aerodynamic force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intialize parameters
global A B Fc Fa

A = Ain;
B = Bin;
Fc = Fcin;
precent = 0;
compare = norm(States);

dxStates = States(1:6*SegNum);
xStates = States(6*SegNum+1:end);
Theta = [xStates(5:6:3*SegNum);0;xStates(3*SegNum+5:6:end)];
PhiTheta = [xStates(6:6:3*SegNum);0;xStates(3*SegNum++6:6:end)];
Z = [xStates(1:6:3*SegNum);0;xStates(3*SegNum+1:6:end)];
Phiz = [xStates(2:6:3*SegNum);0;xStates(3*SegNum+2:6:end)];
X = [xStates(3:6:3*SegNum);0;xStates(3*SegNum+3:6:end)];
Phix = [xStates(4:6:3*SegNum);0;xStates(3*SegNum+4:6:end)];

Fa = zeros(6*SegNum,1);
ref = setRef(L,bw,cord,bc,centroid,0);

while precent < convergePercent
    lattice = generateLattice5_30(SegNum,Z,Phiz,X,Phix,Theta,PhiTheta,cord,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,state);
    geo = setupGeo5_7_2015(zeros(1,3),zeros(1,3),2*max(max(max(lattice.XYZ))),cordNum,SegNum,numSpanB);
    results = solver5_7_2015(state,geo,lattice);
    results = coeff_create5_7_2015(results,lattice,state,ref,geo);
    
    %Determine AeroForces
    for i = 1:SegNum
        %1st pair is Z displacement, 
        %2nd pair is x displacement,
        %3rd pair twist
        Fa(6*(i-1)+2) = sum(results.M(1:SegNum:SegNum*cordNum,1));
        Fa(6*(i-1)+5) = sum(results.M(1:SegNum:SegNum*cordNum,2));
        Fa(6*(i-1)+4) = sum(results.M(1:SegNum:SegNum*cordNum,3));
        Fa(6*(i-1)+1) = sum(results.F(1:SegNum:SegNum*cordNum,3));
        Fa(6*(i-1)+3) = sum(results.F(1:SegNum:SegNum*cordNum,1));
    end
    
    [T,XInt] = ode45(@dynamics,0:dt:dt,States);
    
    States = XInt(2,:)';
    dxStates = States(1:6*SegNum);
    xStates = States(6*SegNum+1:end);
    precent = compare/norm(XInt(2,:));
    compare = norm(XInt(2,:));
end
DT = T(2);
end

function dx = dynamics(t,States)

global A B Fc Fa

dx = A*States+B*(Fa+Fc);
end