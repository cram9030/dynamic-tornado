function [masterState,masterLattice,masterGeo,masterResults,ref,forceTwist] = aeroCoefficient(M,C,K,CG,ActLoc,alphaRange,twistRange,aileronRange,betaRange,airSpeedRange,airDensity,Wings,S_ref,C_ref,B_ref,mac_pos,step)

%alphaRange, twistRange, and aileronRange must be in degrees
%Might only work with the same panel numbers on all wings will need to
%check this later

initTime = cputime;

K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

forceTwist = zeros(length(twistRange),length(alphaRange));

[ref]=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos);
totalPanels = 0;
for h = 1:length(Wings)
    totalPanels = Wings(h).wing.SegNum*Wings(h).wing.cordNum+totalPanels;
end

[numAileron,~] = size(aileronRange);
[numTwist,~] = size(twistRange);

for n = 1:length(airSpeedRange)
    for m = 1:length(betaRange)
        for g = 1:numAileron
            for j = 1:numTwist
                for i = 1:length(alphaRange)
                    
                    [Theta,~,~,~,~,force] = generateStates(K,Wings,[twistRange(j,1)*pi/180,twistRange(j,2)*pi/180]);
                    forceTwist(j,i) = force;
                    
                    masterState(j,i,g,m,n) = setupState2015_8_10(alphaRange(i)*pi/180*ones(totalPanels,1),alphaRange(i)*pi/180,betaRange(m),airSpeedRange(n)*ones(totalPanels,1),airSpeedRange(n),airDensity);
                    masterGeo(j,i,g,m,n) = setupGeo2015_8_10([ActLoc(1),0,ActLoc(2)],[CG(1),0,CG(2)],Wings);
                    for k = 1:length(Wings(1).wing.Controls)
                        masterGeo(j,i,g,m,n).Wings(1).wing.Controls(k).rotation = aileronRange(g,k);
                    end
                    masterGeo(j,i,g,m,n).Wings(1).wing.Theta = Theta;
                    masterGeo(j,i,g,m,n).Wings(1).wing.Flex = [twistRange(j,1)*pi/180,twistRange(j,2)*pi/180];
                    masterLattice(j,i,g,m,n) = generateLattice2015_10_6(masterGeo(j,i,g,m,n),S_ref,C_ref,B_ref,mac_pos,masterState(j,i,g,m,n));
                    
                    geo = masterGeo(j,i,g,m,n);
                    for k = 1:length(masterGeo(j,i,g,m,n).Wings(1).wing.Flex)
                        twist = [twistRange(j,1)*pi/180,twistRange(j,2)*pi/180];
                        twist(k) = twist(k)+1j*step;
                        [Thetai,~,~,~,~,~] = generateStates(K,Wings,twist);
                        geo.Wings(1).wing.Theta = Thetai;
                        latticei(k) = generateLattice2015_10_6(geo,S_ref,C_ref,B_ref,mac_pos,masterState(j,i,g,m,n));
                    end
                    geo.Wings(1).wing.Theta = Theta;
                    for k = 1:length(Wings(1).wing.Controls)
                        geo.Wings(1).wing.Controls(k).rotation = aileronRange(g,k)+1j*step;
                        latticei(k+length(masterGeo(j,i,g,m,n).Wings(1).wing.Flex)) = generateLattice2015_10_6(geo,S_ref,C_ref,B_ref,mac_pos,masterState(j,i,g,m,n));
                    end

                    results = dynamicSolver(masterState(j,i,g,m,n),masterGeo(j,i,g,m,n),masterLattice(j,i,g,m,n),latticei,step);
                    masterResults(j,i,g,m,n) = coeff_create(results,masterLattice(j,i,g,m,n),masterState(j,i,g,m,n),ref,masterGeo(j,i,g,m,n));
                end
                disp(['Precent Complete: ',num2str(100*(length(betaRange)*length(aileronRange)*length(twistRange)*(n-1)+length(aileronRange)*length(twistRange)*(m-1)+length(twistRange)*(g-1)+j)/(length(twistRange)*length(aileronRange)*length(betaRange)*length(airSpeedRange)))]);
                disp(['Elapse Time: ' num2str(cputime-initTime)]);
            end
        end
    end
end
end

function [Theta,Phiz,Phix,X,Z,force] = generateStates(K,Wings,twist)
    States = zeros(5*(Wings(1).wing.SegNum+1),1);
    tempK = zeros(size(K));
    tempTwist = zeros(length(tempK),1);

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
    tempTwist(end) = -K(end,end)*twist(1);
    tempTwist(end-5) = -K(end-5,end)*twist(1);
    staticTwist = tempK\tempTwist;
    staticTwist(end) = twist(1);

    %Convert to passable twist states
    for k = 5:5:5/2*Wings(1).wing.SegNum
        States(k) = staticTwist(length(staticTwist)-k+5);
    end

    tempK(end,end) = -1;
    tempK(end-5,end) = 0;
    tempTwist(end) = -K(end,end)*twist(2);
    tempTwist(end-5) = -K(end-5,end)*twist(2);
    staticTwist = tempK\tempTwist;
    force = staticTwist(end);
    staticTwist(end) = twist(2);

    States(length(staticTwist)+6:end) = staticTwist;
    Theta = States(5:5:5*(Wings(1).wing.SegNum+1));
    Phiz = States(2:5:5*(Wings(1).wing.SegNum+1));
    Phix = States(4:5:5*(Wings(1).wing.SegNum+1));
    X = States(3:5:5*(Wings(1).wing.SegNum+1));
    Z = States(1:5:5*(Wings(1).wing.SegNum+1));
end