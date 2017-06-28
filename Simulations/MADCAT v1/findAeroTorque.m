function aeroTorque = findAeroTorque(Lattice,Results,torqueLoc)

aeroTorque = zeros(length(Results),3);
rotationMatrix = [cos(45*pi/180),-sin(45*pi/180),0;sin(45*pi/180),cos(45*pi/180),0;0,0,1];
for j = 1:length(Results)
    for i = 1:length(Lattice(j).COLLOC)
        if Lattice(j).COLLOC(i,2)>.156/2
            %rotatedLoc = rotationMatrix*(Lattice(j).COLLOC(i,:)-torqueLoc(i,:))';
            rotatedLoc = rotationMatrix*(Lattice(j).COLLOC(i,:)-torqueLoc)';
            rotatedForce = rotationMatrix*Results(j).F(i,:)';
            aeroTorque(j,:) = cross(rotatedLoc',rotatedForce')+aeroTorque(j,:);
        end
    end
end