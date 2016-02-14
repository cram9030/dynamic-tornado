%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,G] = dynamicTwist2016_2_6(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamicTwist2016_2_6: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates the dynamic tip twist of the wing for a
% proscribed tip profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
%           seitupState5_7_2015
%           dynamicSolver6_23_2015
%           setupGeo6_15_2013
%           setRef6_15
% Inputs:   M - wing mass matrix for a single wing side
%           C - wing structural damping matrix for a single wing side
%           K - wing stiffness matrix for a sinlge wing side
%           SegNum - number of wing segments for both wings
%           halfWingLength - half span length
%           cord - wing cord *NOTE: this will need to be adjusted for
%           variable cord wings
%           centroid - centroid of the wing *See note above
%           ActLoc - Actuator location
%           cordNum - number of vortex horse shoes in the cord direction
%           beta - degree of sideslip
%           airSpeed - airspeed in m/2
%           airDensity - Air density in kg/m^3
%           twist - matrix where the first column is a time array and the
%           second is the corisponding tip twist at the specified time
%           Cdp - matrix where the first column is the amount of tiwst and
%           the second column is the parasitic drag coefficent
% Output:   CL - array lift coefficent at the times from twist input
%           CD - array drag coefficent at the times from twist input
%           LD - array lift/drag ratio at the times from twist input
%           forceTwist - array of required torque for specified tip twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/50;
twist(:,1) = x(1:end-1);
twist(:,2) = diff(x)/dt;

%global tunnelWing K M C Cdp S_ref C_ref B_ref alpha_0 centroid ActLoc rhom_inf airSpeed
tunnelWing = [];
airSpeed = [];
cordNum = [];
ActLoc = [];
C_ref = [];
rhom_inf = [];
centroid = [];
S_ref = [];
B_ref = [];

load('OptimizationLoadFile2016_2_8.mat')

%global tunnelWing K M C Cdp S_ref C_ref B_ref alpha_root centroid ActLoc rhom_inf airSpeed

%Intialize clamped boundary condition for structure arrays
K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

%Initalizing output variables
CL = zeros(length(twist),1);
CD = zeros(length(twist),1);

%Intialize length, airspeed, and refrence
ref=setRef2015_8_17(S_ref,C_ref,B_ref,[0,0,0]);
totalPanels = 0;
for h = 1:length(tunnelWing)
    totalPanels = tunnelWing(h).wing.SegNum*tunnelWing(h).wing.cordNum+totalPanels;
end

airfoilRot = zeros(2,2,3);
for i = 1:tunnelWing(1).wing.cordNum
    phi = -atan(ActLoc(2)/(ActLoc(1)-((i-1)*C_ref/tunnelWing(1).wing.cordNum)+.75*C_ref/tunnelWing(1).wing.cordNum));
    airfoilRot(:,:,i) = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end


%Convert alpha_root from degrees to radian
alpha_root = alpha_root*pi/180;

parfor count = 1:length(twist(:,1))
    
    tempWing = tunnelWing;
    %Initialize alphas, Dihedrial angle, and forces to zero
    alpha = zeros(tempWing(h).wing.SegNum,1);
    vLoc = zeros(2,tempWing(h).wing.SegNum);
    alpha_aero = zeros(tempWing(h).wing.SegNum,1);
    artU_infMag = zeros(tempWing(h).wing.SegNum,1);
    
    %Intialize temporary variables
    tempK = zeros(size(K));
    tempTwist = zeros(length(tempK),1);
    tempRotVel = zeros(length(tempK),1);
    States = zeros(5*(tempWing(h).wing.SegNum+1),1);
    
    %Determine twist of wing based off of perscribed tip twist
    tipTwist = twist(count,1);
    tipRotVel = twist(count,2);
    
    tempK = K;
    tempK(end,end) = -1;
    tempK(end-5,end) = 0;
    tempTwist(end) = -K(end,end)*tipTwist;
    tempTwist(end-5) = -K(end-5,end)*tipTwist;
    staticTwist = tempK\tempTwist;
    forceTwist = staticTwist(end);
    staticTwist(end) = tipTwist;
    
    %Determine twist rate from tip twist rate
    tempC = C;
    tempC(end:end,end) = -1;
    tempC(end-5,end) = 0;
    tempRotVel(end) = -C(end,end)*tipRotVel;
    tempRotVel(end-5) = -C(end-5,end)*tipRotVel;
    staticVel = tempC\tempRotVel;
    staticVel(end) = tipRotVel;
    
    %Convert to passable twist states
    for k = 5:5:5/2*tempWing(h).wing.SegNum
        States(k) = staticTwist(length(staticTwist)-k+5);
    end
    States(length(staticTwist)+6:end) = staticTwist;
    Theta = States(5:5:5*(tempWing(h).wing.SegNum+1));
    Phiz = States(2:5:5*(tempWing(h).wing.SegNum+1));
    Phix = States(4:5:5*(tempWing(h).wing.SegNum+1));
    X = States(3:5:5*(tempWing(h).wing.SegNum+1));
    Z = States(1:5:5*(tempWing(h).wing.SegNum+1));
    
    %Calculate local anlge of attack for segments
    for j = 1:tempWing(h).wing.SegNum
        alpha(j) = 0.5*Theta(j+1)+0.5*Theta(j);
    end
    
    %Calculate artificial airspeed and aeroelastic angle of attack
    tempU_inf = [airSpeed*ones(1,tempWing(h).wing.cordNum*(tempWing(h).wing.SegNum));zeros(1,tempWing(h).wing.cordNum*(tempWing(h).wing.SegNum))];
    tempArtU_inf = zeros(2,tempWing(h).wing.SegNum*tempWing(h).wing.cordNum);
    for i = 0:cordNum-1
        tempArtU_inf(:,tempWing(h).wing.SegNum*i+1:tempWing(h).wing.SegNum*(i+1)) = airfoilRot(:,:,i+1)*[zeros(1,tempWing(h).wing.SegNum);norm([ActLoc(1);ActLoc(2)-(i*C_ref/cordNum)+.75*C_ref/cordNum],2)*[flipud(staticVel(5:5:end));staticVel(5:5:end)]'];
    end
    artU_inf = tempArtU_inf;
    alpha = alpha_root*ones(tempWing(h).wing.cordNum*tempWing(h).wing.SegNum,1);
    for i = 1:length(artU_inf)
        vLoc(:,i) = [cos(alpha(i)) -sin(alpha(i));sin(alpha(i)) cos(alpha(i))]*tempU_inf(:,i)+artU_inf(:,i);
        alpha_aero(i) = atan(vLoc(2,i)/vLoc(1,i));
        artU_infMag(i) = norm(vLoc(:,i),2);
    end
    
    tempWing(1).wing.Theta = Theta;
    tempWing(1).wing.Phiz = Phiz;
    tempWing(1).wing.Phix = Phix;
    tempWing(1).wing.X = X;
    tempWing(1).wing.Z = Z;
    tempWing(1).wing.Flex = [twist(count,1),twist(count,1)];
    
    state = setupState2015_8_10(alpha_aero,alpha_root,0,artU_infMag,airSpeed,rhom_inf);
    geo = setupGeo2015_8_10([ActLoc(1),0,ActLoc(2)],[centroid(1),0,centroid(2)],tempWing);
    lattice = generateLattice2015_11_18(geo,S_ref,C_ref,B_ref,[centroid(1),0,centroid(2)],state);
    latticei = repmat(lattice,1);
    results = dynamicSolver(state,geo,lattice,latticei,1e-20,0);
    results = coeff_create(results,lattice,state,ref,geo,0);
    CL(count) = results.CL;
    CD(count) = results.CD+interp1(Cdp(:,1),Cdp(:,2),tipTwist);
    %disp(['Precent Complete: ',num2str(count/length(twist)*100)])
end

F = [1000*mean(CD*.5*rhom_inf*U_inf^2*S_ref);
    mean(CL)*.5*rhom_inf*U_inf^2*S_ref;
    max(abs(twist(:,2)))/dt;
    x(1)-x(end)
    CL*.5*rhom_inf*U_inf^2*S_ref;
    ];
<<<<<<< HEAD
<<<<<<< HEAD
F(1)
=======
=======
>>>>>>> origin/master
F(1)/1000

>>>>>>> f67a079be5e83328ebf301c26d36a9cfc816aadb
% Define the derivatives.
G=[];