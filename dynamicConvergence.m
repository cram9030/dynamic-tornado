function [Fa] = dynamicConvergence(A,B,Fc,SegNum,States,convergePercent,cord,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,state)

precent = 0;
G = inv(A);
Theta = States(5:6:6*(SegNum+1));
PhiTheta = States(6:6:6*(SegNum+1));
Z = States(1:6:6*(SegNum+1));
Phiz = States(2:6:6*(SegNum+1));
X = States(3:6:6*(SegNum+1));
Phix = States(4:6:6*(SegNum+1));

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
    
    compare = [G(6*SegNum+1:end,6*SegNum+1:end),G(6*SegNum+1:end,1:6*SegNum)-G(6*SegNum+1:end,6*SegNum+1:end)*G(1:6*SegNum,1:6*SegNum)-eye(6*SegNum)]*States;
    precent = norm(compare)/norm(B*(Fa+Fc));
end