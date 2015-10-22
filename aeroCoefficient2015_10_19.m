function [masterState,masterLattice,masterGeo,masterResults,ref,forceTwist] = aeroCoefficient2015_10_19(M,C,K,CG,ActLoc,alphaRange,twistRange,aileronRange,beta,airSpeed,airDensity,Wings,S_ref,C_ref,B_ref,mac_pos)

%alphaRange, twistRange, and aileronRange must be in degrees
%Might only work with the same panel numbers on all wings will need to
%check this later

K = K(6:end,6:end);
M = M(6:end,6:end);
C = C(6:end,6:end);

forceTwist = zeros(length(twistRange),length(alphaRange));

[ref]=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos);
totalPanels = 0;
for h = 1:length(Wings)
    totalPanels = Wings(h).wing.SegNum*Wings(h).wing.cordNum+totalPanels;
end

for g = 1:length(aileronRange)
    for j = 1:length(twistRange)
        parfor i = 1:length(alphaRange)
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
            tempTwist(end) = -K(end,end)*twistRange(j,1)*pi/180;
            tempTwist(end-5) = -K(end-5,end)*twistRange(j,1)*pi/180;
            staticTwist = tempK\tempTwist;
            forceTwist(j,i) = staticTwist(end);
            staticTwist(end) = twistRange(j,1)*pi/180;

            %Convert to passable twist states
            for k = 5:5:5/2*Wings(1).wing.SegNum
                States(k) = staticTwist(length(staticTwist)-k+5);
            end

            tempK(end,end) = -1;
            tempK(end-5,end) = 0;
            tempTwist(end) = -K(end,end)*twistRange(j,2)*pi/180;
            tempTwist(end-5) = -K(end-5,end)*twistRange(j,2)*pi/180;
            staticTwist = tempK\tempTwist;
            forceTwist(j,i) = staticTwist(end);
            staticTwist(end) = twistRange(j,2)*pi/180;
            
            States(length(staticTwist)+6:end) = staticTwist;
            Theta = States(5:5:5*(Wings(1).wing.SegNum+1));
            Phiz = States(2:5:5*(Wings(1).wing.SegNum+1));
            Phix = States(4:5:5*(Wings(1).wing.SegNum+1));
            X = States(3:5:5*(Wings(1).wing.SegNum+1));
            Z = States(1:5:5*(Wings(1).wing.SegNum+1));
            
%             for k = 1:length(Wings(1).wing.Controls)
%                 Wings(1).wing.Controls(k).rotation = aileronRange(g,k)*pi/180;
%             end

            masterState(j,i,g) = setupState2015_8_10(alphaRange(i)*pi/180*ones(totalPanels,1),alphaRange(i)*pi/180,beta,airSpeed*ones(totalPanels,1),airSpeed,airDensity);
            masterGeo(j,i,g) = setupGeo2015_8_10([ActLoc(1),0,ActLoc(2)],[CG(1),0,CG(2)],Wings);
            for k = 1:length(Wings(1).wing.Controls)
                masterGeo(j,i,g).Wings(1).wing.Controls(k).rotation = aileronRange(g,k);
            end
            masterGeo(j,i,g).Wings(1).wing.Theta = Theta;
            masterLattice(j,i,g) = generateLattice2015_10_6(masterGeo(j,i,g),S_ref,C_ref,B_ref,mac_pos,masterState(j,i,g));
            results = dynamicSolver6_23_2015(masterState(j,i,g),masterGeo(j,i,g),masterLattice(j,i,g));
            masterResults(j,i,g) = coeff_create2015_8_13(results,masterLattice(j,i,g),masterState(j,i,g),ref,masterGeo(j,i,g));
        end
        disp(['Precent Complete: ',num2str(100*(length(twistRange)*j*(g-1)+j)/(length(twistRange)*length(aileronRange)))]);
        disp(['CPU Time: ' num2str(cputime)]);
    end
end