function [polyCoef,spanCoeff] = generatePoly(geo,S_ref,C_ref,B_ref,mac_pos,AoA,speed,airDensity,totalPanels,twist,spanDegree,angleDegree)

%Intialize output
liftPolyCoef = zeros(length(AoA),spanDegree+1);
polyCoef = zeros(spanDegree+1,2*(angleDegree+1));
spanCoeff = zeros(length(twist),angleDegree+1,spanDegree+1);
a = zeros(length(twist),2);
b = zeros(length(twist),2);
c = zeros(length(twist),2);
d = zeros(length(twist),2);
e = zeros(length(twist),2);

%Iterate through all the requested twists
for i = 1:length(twist)
    [spanCL,~,~,~,Results,~,~] = symmetricAoASweep(geo,S_ref,C_ref,B_ref,mac_pos,AoA,speed,airDensity,totalPanels,twist(i));
    % Find coefficient for fourth order polynomial for each spanwise lift
    % coefficents
    for j = 1:length(AoA)
        liftPolyCoef(j,:) = polyfit(Results(j).spanResults.ystation*39.3701,spanCL(:,j),spanDegree);
    end
    % Fit the spanwise coefficent linear fit for the angle of attack
    for j = 1:spanDegree+1
        spanCoeff(i,:,j) = polyfit(AoA',liftPolyCoef(:,j),angleDegree);
    end
end

% Fit the angle of attack coefficents to twist coefficents
for i = 1:spanDegree+1 
    for j = 1:angleDegree+1
        polyCoef(i,(j-1)*(angleDegree+1)+1:j*(angleDegree+1)) = polyfit(twist',spanCoeff(:,j,i),angleDegree);
    end
end