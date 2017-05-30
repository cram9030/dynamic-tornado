%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CL_local = polyCheck(polyCoef,yStation,AoA,twist,spanDegree,angleDegree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyCheck: Function for MADCAT v1 Simulations					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is intended to check the results of the generatePoly
%function that is used to generate a polynomial estimate of the spanwise
%lifting values with respect to the angle of attack and tip twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   polyCoef -  this is a 4X5 Matrix where for each column the pair
%                       of rows are a linear polynomial representation of the angle of attack
%                       coefficents representing the spanwise lifting distribution
% Output:   CL_local -  the resulting local lift coefficents using the
%                       input polynomical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intialize local lift parameters
CL_local = zeros(length(AoA),length(yStation));

%Iterate through the angle of attackes
for i = 1:length(AoA)
    abAlpha = zeros(spanDegree+1,angleDegree+1);
    %Use the polynomial fitting of the twist to get the coefficient of the
    %polynomial fit of the spanwise lift coefficents with angle of attack
    for j = 1:spanDegree+1
        for k = 1:angleDegree+1
            abAlpha(j,k) = polyval(polyCoef(j,(k-1)*(angleDegree+1)+1:k*(angleDegree+1)),twist);
        end
    end
    liftPoly = zeros(1,spanDegree+1);
    %Use the polynomial fitting of angle of attack to the coefficents of
    %the spanwise lifting distribution
    for j = 1:spanDegree+1
        liftPoly(j) = polyval(abAlpha(j,:),AoA(i));
    end
    %Polpulate the spanwise lifting coefficent
    CL_local(i,:) = polyval(liftPoly,yStation);
end
CL_local = CL_local';