%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meanCord] = createMeanCord(filename,cordNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setupState5_7_2015: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the mean cord line for a given airfoil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
% Inputs:   filename - filename containing 
% Output:   meanCord - mean cord array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import airfoil
airfoil = importdata(filename);
try
    data = airfoil.data;
catch
    data = airfoil;
end

%Find split between top and bottom and top of the airfoil
zeroLoc = find(data(:,1) == 0);

%Get top and bottom cord splits then determine mean cord
xCord = 0:1/cordNum:1;
yCordTop = interp1(data(1:zeroLoc,1),data(1:zeroLoc,2),xCord);
yCordBot = interp1(data(zeroLoc:end,1),data(zeroLoc:end,2),xCord);
meanCord = [xCord',(yCordTop'-yCordBot')/2+yCordBot'];