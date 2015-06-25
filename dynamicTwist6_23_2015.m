%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [masterState,masterLattice,masterGeo,masterResults,ref,CL,CD,LD,forceTwist] = dynamicTwist6_23_2015(M,C,K,SegNum,halfWingLength,cord,centroid,ActLoc,cordNum,bw,bc,noseH,zb,alpha_root,numSpanB,beta,airSpeed,airDensity,twist,Cdp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamicTwist6_23_2015: Function for Dynamic TORNADO						
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

%Intialize clamped boundary condition for structure arrays
K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

%Intialize temporary variables
tempK = zeros(size(K));
tempTwist = zeros(length(tempK),1);
tempRotVel = zeros(length(tempK),1);
States = zeros(5*(SegNum+1),1);

%Initalizing output variables
forceTwist = zeros(length(twist),1);
CL = zeros(length(twist),1);
CD = zeros(length(twist),1);
LD = zeros(length(twist),1);

%Intialize length and refrence
L = 2*halfWingLength/SegNum*ones(SegNum,1);
[ref]=setRef6_15(L,cord,zeros(1,3));

%Convert alpha_root from degrees to radian
alpha_root = alpha_root*pi/180;

count = 1;
for t = twist(:,1)'
    %Initialize alphas, Dihedrial angle, and forces to zero
    DiHiAng = zeros(SegNum,1);
    alpha = zeros(SegNum,1);
    
    %Determine twist of wing based off of perscribed tip twist
    tipTwist = twist(count,2);
    tipRotVel = twist(count,3);
    
    
    tempK = K;
    tempK(end,end) = -1;
    tempK(end-5,end) = 0;
    tempTwist(end) = -K(end,end)*tipTwist;
    tempTwist(end-5) = -K(end-5,end)*tipTwist;
    staticTwist = tempK\tempTwist;
    forceTwist(count) = staticTwist(end);
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
    for k = 5:5:5/2*SegNum
        States(k) = staticTwist(length(staticTwist)-k+5);
    end
    States(length(staticTwist)+6:end) = staticTwist;
    Theta = States(5:5:5*(SegNum+1));
    Phiz = States(2:5:5*(SegNum+1));
    Phix = States(4:5:5*(SegNum+1));
    X = States(3:5:5*(SegNum+1));
    Z = States(1:5:5*(SegNum+1));
    
    %Calculate local anlge of attack for segments
    for j = 1:SegNum
        alpha(j) = 0.5*Theta(j+1)+0.5*Theta(j);
    end
    alpha = alpha+alpha_root;
    
    %Calculate artificial airspeed and aeroelastic angle of attack
    tempArtU_inf = sqrt((ActLoc(1)*[flipud(staticVel(5:5:end));staticVel(5:5:end)]+airSpeed*sin(alpha)).^2+(airSpeed*cos(alpha)).^2);
    artU_inf = [tempArtU_inf(1:SegNum/2);airSpeed*ones(numSpanB,1);tempArtU_inf(SegNum/2+1:end)];
    tempAlpha_aero = atan((ActLoc(1)^2*[flipud(staticVel(5:5:end));staticVel(5:5:end)]+airSpeed*sin(alpha))./(airSpeed*cos(alpha)));
    alpha_aero = [tempAlpha_aero(1:SegNum/2);alpha_root*ones(numSpanB,1);tempAlpha_aero(SegNum/2+1:end)];
    
    masterState(count) = setupState6_24_2015(alpha_aero,alpha_root,beta,artU_inf,airSpeed,airDensity);
    masterGeo(count) = setupGeo6_15_2015([ActLoc(1),0,ActLoc(2)],[centroid(1),0,centroid(2)],sum(L),cordNum,SegNum);
    masterLattice(count) = generateLattice6_24_2015(SegNum,Z,Phiz,X,Phix,Theta,cord,masterGeo(count),cordNum,L,bw,bc,noseH,zb,numSpanB,masterState(count));
    results = dynamicSolver6_23_2015(masterState(count),masterGeo(count),masterLattice(count));
    masterResults(count) = coeff_create6_24_2015(results,masterLattice(count),masterState(count),ref,masterGeo(count));
    CL(count) = masterResults(count).CL;
    CD(count) = masterResults(count).CD+interp1(Cdp(:,1),Cdp(:,2),tipTwist);
    LD(count) = CL(count)/CD(count);
    disp(['Precent Complete: ',num2str(count/length(twist)*100)])
    count = count + 1;
end