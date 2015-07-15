function [masterState,masterLattice,masterGeo,masterResults,ref,forceTwist,CL,CD] = staticAlphaSweepBi7_10_2015(M,C,K,SegNum,halfWingLength,cord,centroid,ActLoc,cordNum,bw,bc,bh,noseH,zb,numSpanB,alphaRange,twistRange,beta,airSpeed,airDensity)

K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

tempK = zeros(size(K));
tempTwist = zeros(length(tempK),1);
States = zeros(5*(SegNum+1),1);

forceTwist = zeros(length(twistRange),length(alphaRange));
CL = zeros(length(twistRange),length(alphaRange));
CD = zeros(length(twistRange),length(alphaRange));

%Intialize length, airspeed, and refrence
L = 2*halfWingLength/SegNum*ones(SegNum,1);
[ref]=setRef6_15(L,cord,zeros(1,3));
phi = -atan(ActLoc(1)/ActLoc(2));
airfoilRot = [cos(phi) -sin(phi);sin(phi) cos(phi)];

for j = twistRange
    for i = alphaRange
        %Determine twist of wing based off of perscribed tip twist
        tempK = K;
        %Will be switching the states around to solve for the twist and
        %force of the sections. Instead of F = KX where X is
        %[z_i,phi_zi,x_i,phi_xi,theta_i,...,z_end,phi_zend,x_end,phi_xend,theta_end]
        %it will end with ,z_end,phi_zend,x_end,phi_xend,torque_tip]
        %meaning that the stiffness matrix must have the components that
        %are dependent on the tip twist altered. The F must then have the
        %corrisponding components replaced with the stiffness times the tip
        %twist then inverted.
        tempK(end,end) = -1;
        tempK(end-5,end) = 0;
        tempTwist(end) = -K(end,end)*j*pi/180;
        tempTwist(end-5) = -K(end-5,end)*j*pi/180;
        staticTwist = tempK\tempTwist;
        forceTwist(j-min(twistRange)+1,i-min(alphaRange)+1) = staticTwist(end);
        staticTwist(end) = j*pi/180;
        
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
        alpha = zeros(SegNum,1);
        for k = 1:SegNum
            alpha(k) = 0.5*Theta(k+1)+0.5*Theta(k);
        end
        alpha = i*pi/180*ones(SegNum,1)-alpha;
        alpha = [alpha;i*pi/180*ones(4*numSpanB,1)];
        
        masterState(j-min(twistRange)+1,i-min(alphaRange)+1) = setupState6_24_2015(alpha,i*pi/180,beta,airSpeed*ones(size(alpha)),airSpeed,airDensity);
        masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1) = setupGeo6_15_2015([ActLoc(1),0,ActLoc(2)],[centroid(1),0,centroid(2)],sum(L),cordNum,SegNum);
        masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1) = generateLatticeBi7_10_2015(SegNum,Z,Phiz,X,Phix,Theta,cord,masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1),cordNum,L,bw,bc,bh,noseH,zb,numSpanB,masterState(j-min(twistRange)+1,i-min(alphaRange)+1));
        results = dynamicSolver6_23_2015(masterState(j-min(twistRange)+1,i-min(alphaRange)+1),masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1),masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1));
        masterResults(j-min(twistRange)+1,i-min(alphaRange)+1) = coeff_create6_24_2015(results,masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1),masterState(j-min(twistRange)+1,i-min(alphaRange)+1),ref,masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1));
        CL(j-min(twistRange)+1,i-min(alphaRange)+1) = masterResults(j-min(twistRange)+1,i-min(alphaRange)+1).CL;
        CD(j-min(twistRange)+1,i-min(alphaRange)+1) = masterResults(j-min(twistRange)+1,i-min(alphaRange)+1).CD;
        disp(['Precent Complete: ',num2str((length(alphaRange)*(j-min(twistRange))+i-min(alphaRange)+1)/(length(twistRange)*length(alphaRange)))])
    end
end