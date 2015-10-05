%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = createCPMovie(lattice,results,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createCPMovie: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates the dynamic tip twist of the wing for a
% proscribed tip profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	MATLAB 5.2 std fcns
%
% Inputs:   lattice - lattice structure
%           results - array of result structures
%           filename - name of video file that will be saved
% Output:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(results);

F(N) = struct('cdata',[],'colormap',[]);

figure(4)% 	Delta cp plot
rotate3d on 
colormap(hot);
fill3(lattice(1).XYZ(:,:,1)',lattice(1).XYZ(:,:,2)',lattice(1).XYZ(:,:,3)',results(1).cp')
minLattice = min(min([lattice(:).XYZ]));
maxLattice = max(max([lattice(:).XYZ]));
axis equal
axis([minLattice(1) maxLattice(1) minLattice(2) maxLattice(2) minLattice(3) maxLattice(3)])
caxis([min(min([results(:).cp]')) max(max([results(:).cp]'))])
title('C_p Distribution')
colorbar('vert')
ax = gca;
ax.NextPlot = 'replaceChildren';

for j = 1:N
    fill3(lattice(j).XYZ(:,:,1)',lattice(j).XYZ(:,:,2)',lattice(j).XYZ(:,:,3)',results(j).cp')
    drawnow
    F(j) = getframe(gcf);
end

movie2avi(F,filename,'compression','None')