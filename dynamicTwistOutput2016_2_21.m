%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CL,CD,Cm,Cl,Cn,LD] = dynamicTwistOutput2016_2_21(x)
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
% Inputs:  
% Output:   CL - array lift coefficent at the times from twist input
%           CD - array drag coefficent at the times from twist input
%           LD - array lift/drag ratio at the times from twist input
%           forceTwist - array of required torque for specified tip twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tunnelWing = [];
airSpeed = [];
cordNum = [];
ActLoc = [];
C_ref = [];
rhom_inf = [];
centroid = [];
S_ref = [];
B_ref = [];
alpha_root = [];
sideSlip = [];

load('OptimizationLoadFile2016_2_21.mat')

dt = 1/50;
twist(:,1) = x(1:length(x)/2);
twist(:,2) = x(1+length(x)/2:end);
twistRate(:,1) = diff(twist(:,1))/dt;
twistRate(:,2) = diff(twist(:,2))/dt;

CL = zeros(length(twist(:,1))-1,1);
CD = zeros(length(twist(:,1))-1,1);
Cm = zeros(length(twist(:,1))-1,1);
Cl = zeros(length(twist(:,1))-1,1);
Cn = zeros(length(twist(:,1))-1,1);

%Intialize clamped boundary condition for structure arrays
K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

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

parfor count = 1:length(twist(:,1))-1
    
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
    
    %Determine Left twist of wing based off of perscribed tip twist
    tipTwistL = twist(count,1);
    tipRotVelL = twistRate(count,1);
    
    tempK = K;
    tempK(end,end) = -1;
    tempK(end-5,end) = 0;
    tempTwist(end) = -K(end,end)*tipTwistL;
    tempTwist(end-5) = -K(end-5,end)*tipTwistL;
    staticTwistL = tempK\tempTwist;
    staticTwistL(end) = tipTwistL;
    
    %Determine left twist rate from tip twist rate
    tempC = C;
    tempC(end:end,end) = -1;
    tempC(end-5,end) = 0;
    tempRotVel(end) = -C(end,end)*tipRotVelL;
    tempRotVel(end-5) = -C(end-5,end)*tipRotVelL;
    staticVelL = tempC\tempRotVel;
    staticVelL(end) = tipRotVelL;
    
    %Determine right twist of wing based off of perscribed tip twist
    tipTwistR = twist(count,2);
    tipRotVelR = twistRate(count,2);
    
    tempK = K;
    tempK(end,end) = -1;
    tempK(end-5,end) = 0;
    tempTwist(end) = -K(end,end)*tipTwistR;
    tempTwist(end-5) = -K(end-5,end)*tipTwistR;
    staticTwistR = tempK\tempTwist;
    staticTwistR(end) = tipTwistR;
    
    %Determine right twist rate from tip twist rate
    tempC = C;
    tempC(end:end,end) = -1;
    tempC(end-5,end) = 0;
    tempRotVel(end) = -C(end,end)*tipRotVelR;
    tempRotVel(end-5) = -C(end-5,end)*tipRotVelR;
    staticVelR = tempC\tempRotVel;
    staticVelR(end) = tipRotVelR;
    
    %Convert to passable twist states
    for k = 5:5:5/2*tempWing(h).wing.SegNum
        States(k) = staticTwistL(length(staticTwistL)-k+5);
    end
    States(length(staticTwistR)+6:end) = staticTwistR;
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
        tempArtU_inf(:,tempWing(h).wing.SegNum*i+1:tempWing(h).wing.SegNum*(i+1)) = airfoilRot(:,:,i+1)*[zeros(1,tempWing(h).wing.SegNum);norm([ActLoc(1);ActLoc(2)-(i*C_ref/cordNum)+.75*C_ref/cordNum],2)*[flipud(staticVelL(5:5:end));staticVelR(5:5:end)]'];
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
    
    state = setupState2015_8_10(alpha_aero,alpha_root,sideSlip,artU_infMag,airSpeed,rhom_inf);
    geo = setupGeo2015_8_10([ActLoc(1),0,ActLoc(2)],[centroid(1),0,centroid(2)],tempWing);
    lattice = generateLattice2015_11_18(geo,S_ref,C_ref,B_ref,[centroid(1),0,centroid(2)],state);
    latticei = repmat(lattice,1);
    results = dynamicSolver(state,geo,lattice,latticei,1e-20,0);
    results = coeff_create(results,lattice,state,ref,geo,0);
    CL(count) = results.CL;
    CD(count) = results.CD+interp1(Cdp(:,1),Cdp(:,2),tipTwistR);
    Cm(count) = results.Cm;
    Cl(count) = results.Cl;
    Cn(count) = results.Cn;
    %disp(['Precent Complete: ',num2str(count/length(twist)*100)])
end

LD = CL./CD;