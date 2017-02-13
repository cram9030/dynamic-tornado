function [masterState,masterLattice,masterGeo,masterResults,ref,forceTwist,CL,CD] = staticAlphaSweep(M,C,K,SegNum,halfWingLength,cord,centroid,ActLoc,cordNum,alphaRange,twistRange,beta,airSpeed,airDensity)

K = K(7:end,7:end);
M = M(7:end,7:end);
C = C(7:end,7:end);

tempK = zeros(size(K));
tempTwist = zeros(length(tempK),1);
States = zeros(6*(SegNum+1),1);

forceTwist = zeros(length(twistRange),length(alphaRange));
CL = zeros(length(twistRange),length(alphaRange));
CD = zeros(length(twistRange),length(alphaRange));

L = 2*halfWingLength/SegNum*ones(SegNum,1);
[ref]=setRef6_15(L,cord,zeros(1,3));

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
        tempK(end-1:end,end-1) = [-1;0];
        tempK(end-7:end-6,end-1) = [0;0];
        tempTwist(end-1:end) = -K(end-1:end,end-1)*j*pi/180;
        tempTwist(end-7:end-6) = -K(end-7:end-6,end-1)*j*pi/180;
        staticTwist = tempK\tempTwist;
        forceTwist(j-min(twistRange)+1,i-min(alphaRange)+1) = staticTwist(end-1);
        staticTwist(end-1) = j*pi/180;
        
        %Convert to passable twist states
        for k = 5:6:3*SegNum
            States(k:k+1) = staticTwist(length(staticTwist)-1-k+5:length(staticTwist)-k+5);
        end
        States(length(staticTwist)+7:end) = staticTwist;
        Theta = States(5:6:6*(SegNum+1));
        PhiTheta = States(6:6:6*(SegNum+1));
        Phiz = States(2:6:6*(SegNum+1));
        Phix = States(4:6:6*(SegNum+1));
        X = States(3:6:6*(SegNum+1));
        Z = States(1:6:6*(SegNum+1));
        
        masterState(j-min(twistRange)+1,i-min(alphaRange)+1) = setupState5_7_2015(i,beta,airSpeed,airDensity);
        masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1) = setupGeo6_15_2015([ActLoc(1),0,ActLoc(2)],[centroid(1),0,centroid(2)],sum(L),cordNum,SegNum);
        masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1) = generateLatticeNoBody6_15(SegNum,Z,Phiz,X,Phix,Theta,PhiTheta,cord,masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1),cordNum,L,masterState(j-min(twistRange)+1,i-min(alphaRange)+1));
        results = solver5_7_2015(masterState(j-min(twistRange)+1,i-min(alphaRange)+1),masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1),masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1));
        masterResults(j-min(twistRange)+1,i-min(alphaRange)+1) = coeff_create6_15_2015(results,masterLattice(j-min(twistRange)+1,i-min(alphaRange)+1),masterState(j-min(twistRange)+1,i-min(alphaRange)+1),ref,masterGeo(j-min(twistRange)+1,i-min(alphaRange)+1));
        CL(j-min(twistRange)+1,i-min(alphaRange)+1) = masterResults(j-min(twistRange)+1,i-min(alphaRange)+1).CL;
        CD(j-min(twistRange)+1,i-min(alphaRange)+1) = masterResults(j-min(twistRange)+1,i-min(alphaRange)+1).CD;
        disp(['Precent Complete: ',num2str((length(alphaRange)*(j-min(twistRange))+i-min(alphaRange)+1)/(length(twistRange)*length(alphaRange)))])
    end
end