function [results]=dynamicSolver_analytic(state,geo,lattice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 1999, 2007 Tomas Melin
%
% This file is part of Tornado
%
% Tornado is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2, or (at your option) any later version.
%
% Tornado is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public
% License along with Tornado; see the file GNU GENERAL 
% PUBLIC LICENSE.TXT.  If not, write to the Free Software 
% Foundation, 59 Temple Place -Suite 330, Boston, MA
% 02111-1307, USA.
%
% usage: [RESULTS] = solver8(results,state,geo,lattice,ref)
%
% This function computes forces and moments on each panel.
% Inputs are coordinades for old resluts, collocationpoints, 
%   vorticies and Normals, reference area and chord
%
% Example:
%
%   [results]=solver8(results,state,geo,lattice,ref);
%e 
% Calls:
%           Setboundary    
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado core function
%
% Revision History:
%   Bristol,  2007 06 27:  Addition of new header. TM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intializing Results
results = [];

[a vor_length void]=size(lattice.VORTEX);%extracting number of sections in 
										   %"horseshoes"

[w2 void]=fastdw(lattice);

results.dwcond=cond(w2);
%disp('dnwash... ok')
%count=flops         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up right hand side %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs = setDynamicBoundary6_23_2015(lattice,state,geo)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving for rhs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

gamma=w2\rhs';
%disp('gauss... ok')




if state.pgcorr==1
    %tdisp('Trying PG correction')
    %Prandtl-Glauert compressibility correction
    [state.rho sos p_1]=ISAtmosphere(state.ALT);
    M=state.AS/sos;
    corr=1/(sqrt(1-M^2));  
    gamma=gamma*corr;
    %Yeah, this part is not validated yet... or even published, but it
    %seems to work. Do use it with caution though as the math rigour
    %isnt there yet.   
end




b1=vor_length/2;

p1(:,:)=lattice.VORTEX(:,b1,:);		%Calculating panel vortex midpoint	
p2(:,:)=lattice.VORTEX(:,b1+1,:);	%to use as a force locus
lattice.COLLOC(:,:)=(p1+p2)./2;	    % LOCAL control point, vortex midpoint.

c3=lattice.COLLOC-ones(size(lattice.COLLOC,1),1)*geo.CG;

[w3 DW]=fastdw(lattice);	                    %Calculating downwash on vorticies
w4=sum(DW,2);					                %superpositioning aerodynamic influence

DWX=DW(:,:,1);
DWY=DW(:,:,2);
DWZ=DW(:,:,3);

[void nofderiv]=size(gamma);
le=(p2-p1);					%Vortex span vector

Lle=sqrt(sum(le.^2,2));
lehat(:,1)=le(:,1)./Lle;
lehat(:,2)=le(:,2)./Lle;
lehat(:,3)=le(:,3)./Lle;


for j=1:nofderiv
    IW(:,j,1)=DWX*gamma(:,j);
    IW(:,j,2)=DWY*gamma(:,j);
    IW(:,j,3)=DWZ*gamma(:,j);
    
    G(:,1)=gamma(:,j).*lehat(:,1);	%Aligning vorticity along panel vortex
    G(:,2)=gamma(:,j).*lehat(:,2);
    G(:,3)=gamma(:,j).*lehat(:,3);

    V=state.AS;
    Wind=([V.*cos(state.alpha).*cos(state.betha) -V.*cos(state.alpha).*sin(state.betha) V.*sin(state.alpha)])-reshape(IW(:,j,:),[length(IW),3]);
    for i=1:a
        %Wind(i,:)=wind-squeeze(IW(i,j,:))';
        Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]); %Calculating rotations
    end                                   %^^^^^^^---new stuff in T131         %Thanks Luca for pointing out the error here

    Wind=Wind+Rot;								%Adding rotations
    Fprim(:,j,:)=state.rho*cross(Wind,G);			    %Force per unit length
    

    F(:,j,1)=Fprim(:,j,1).*Lle;				%Force per panel
    F(:,j,2)=Fprim(:,j,2).*Lle;				%Force per panel
    F(:,j,3)=Fprim(:,j,3).*Lle;				%Force per panel
    

    C3(:,:,1)=c3(:,1)*ones(1,nofderiv);
    C3(:,:,2)=c3(:,2)*ones(1,nofderiv); 
    C3(:,:,3)=c3(:,3)*ones(1,nofderiv); 
        
end

invW2 = inv(w2);
results.dF_a = zeros(length(F),3,length(F));
results.dF_b = zeros(length(F),3,length(F));
for i = 1:length(F)

    results.dF_a(:,:,i) = state.AS(i)^2*state.rho*[invW2(:,i).*(-(lehat(i,3)*(cos(state.alpha(i))*sin(state.betha) + DWY*lehat(i,2)*invW2(:,i).*...
    (lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha))...
    - lehat(i,3)*(sin(state.alpha(i))*sin(state.betha) - DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3) *cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i))...
    + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha))...
    + lehat(i,2)*(cos(state.alpha(i)) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))...
    *(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha))...
    + lehat(i,2)*(sin(state.alpha(i)) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))),...
    ...
    invW2(:,i).*((lehat(i,3)*(cos(state.betha)*sin(state.alpha(i)) + DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))...
    *(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)) - lehat(i,3)*(cos(state.alpha(i))*cos(state.betha)...
    - DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha(i))...
    - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)) + lehat(i,1)*(cos(state.alpha(i)) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha(i))...
    - lehat(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)...
    - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)) + lattice.N(i,1)*(sin(state.alpha(i)) + DWZ*lattice.N(i,3)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)...
    - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))),...
    ...
    invW2(:,i).*((lehat(i,2)*(cos(state.alpha(i))*cos(state.betha) - DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)) - lehat(i,2)*(cos(state.betha)*sin(state.alpha(i))...
    + DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i)) + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha(i))...
    + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)) + lehat(i,1)*(cos(state.alpha(i))*sin(state.betha) + DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha(i))...
    + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i))...
    + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)) - lattice.N(i,1)*(sin(state.alpha(i))*sin(state.betha) - DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha(i)) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i))...
    + lattice.N(i,2)*sin(state.alpha(i))*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha(i)) + lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha) - lattice.N(i,2)*cos(state.alpha(i))*sin(state.betha))))];

%     results.dF_a(:,:,i) = invW2(:,i)*state.AS(i)^2*state.rho*[-(lehat(i,1)*lattice.N(i,3)*sin(2*state.alpha(i))+lehat(i,2)*lattice.N(i,2)*sin(state.alpha(i))^2*sin(state.betha)- lehat(i,3)*lattice.N(i,3)*sin(state.alpha(i))^2*sin(state.betha) + lehat(i,2)*lattice.N(i,1)*cos(state.alpha(i))^2*cos(state.betha) - lehat(i,2)*lattice.N(i,1)*cos(state.betha)*sin(state.alpha(i))^2 - lehat(i,2)*lattice.N(i,2)*cos(state.alpha(i))^2*sin(state.betha) + lehat(i,3)*lattice.N(i,3)*cos(state.alpha(i))^2*sin(state.betha) + 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha(i))*sin(state.alpha(i))*sin(state.betha)^2 - 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)*sin(state.alpha(i))*sin(state.betha)),...
%         lehat(i,1)*lattice.N(i,3)*sin(2*state.alpha(i)) - lehat(i,1)*lattice.N(i,1)*cos(state.betha) + lehat(i,3)*lattice.N(i,3)*cos(state.betha) + lehat(i,1)*lattice.N(i,2)*sin(state.betha) + 2*lehat(i,1)*lattice.N(i,1)*cos(state.alpha(i))^2*cos(state.betha) - 2*lehat(i,3)*lattice.N(i,3)*cos(state.alpha(i))^2*cos(state.betha) - 2*lehat(i,1)*lattice.N(i,2)*cos(state.alpha(i))^2*sin(state.betha) + 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)^2*sin(state.alpha(i)) - 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha(i))*cos(state.betha)*sin(state.alpha(i))*sin(state.betha),...
%         -((lehat(i,2)*cos(state.betha) + lehat(i,1)*sin(state.betha))*(lattice.N(i,3) - 2*lattice.N(i,3)*cos(state.alpha(i))^2 + 2*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)*sin(state.alpha(i)) - 2*lattice.N(i,2)*cos(state.alpha(i))*sin(state.alpha(i))*sin(state.betha)))];
    results.dF_b(:,:,i) = invW2(:,i)*state.AS(i)^2*state.rho*[cos(state.alpha(i))*(lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i)) + lehat(i,2)*lattice.N(i,2)*cos(state.betha)*sin(state.alpha(i)) - lehat(i,3)*lattice.N(i,3)*cos(state.betha)*sin(state.alpha(i)) + lehat(i,2)*lattice.N(i,1)*sin(state.alpha(i))*sin(state.betha) - 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)^2 + 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha(i))*cos(state.betha)*sin(state.betha)),...
        cos(state.alpha(i))*(lehat(i,3)*lattice.N(i,2)*cos(state.alpha(i)) + lehat(i,1)*lattice.N(i,2)*cos(state.betha)*sin(state.alpha(i)) + lehat(i,1)*lattice.N(i,1)*sin(state.alpha(i))*sin(state.betha) - lehat(i,3)*lattice.N(i,3)*sin(state.alpha(i))*sin(state.betha) - 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha(i))*cos(state.betha)^2 - 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)*sin(state.betha)),...
        cos(state.alpha(i))*(lehat(i,3)*lattice.N(i,1)*cos(state.alpha(i)) - lehat(i,2)*lattice.N(i,2)*cos(state.alpha(i)) - lehat(i,1)*lattice.N(i,3)*cos(state.betha)*sin(state.alpha(i)) + lehat(i,2)*lattice.N(i,3)*sin(state.alpha(i))*sin(state.betha) - 2*lehat(i,1)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)^2 + 2*lehat(i,1)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)^2 + 2*lehat(i,2)*lattice.N(i,1)*cos(state.alpha(i))*cos(state.betha)*sin(state.betha) + 2*lehat(i,1)*lattice.N(i,2)*cos(state.alpha(i))*cos(state.betha)*sin(state.betha))];
    results.dM_a(:,:,i)=cross(squeeze(C3(:,1,:)),results.dF_a(:,:,i));
    results.dM_b(:,:,i)=cross(squeeze(C3(:,1,:)),results.dF_b(:,:,i));
end

% results.dF_ar = ones(length(w2),1)*[-(lehat(i,1)*lattice.N(i,3)*sin(2*state.alpha_root)+lehat(i,2)*lattice.N(i,2)*sin(state.alpha_root)^2*sin(state.betha)- lehat(i,3)*lattice.N(i,3)*sin(state.alpha_root)^2*sin(state.betha) + lehat(i,2)*lattice.N(i,1)*cos(state.alpha_root)^2*cos(state.betha) - lehat(i,2)*lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root)^2 - lehat(i,2)*lattice.N(i,2)*cos(state.alpha_root)^2*sin(state.betha) + lehat(i,3)*lattice.N(i,3)*cos(state.alpha_root)^2*sin(state.betha) + 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha_root)*sin(state.alpha_root)*sin(state.betha)^2 - 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha)*sin(state.alpha_root)*sin(state.betha)),...
%         lehat(i,1)*lattice.N(i,3)*sin(2*state.alpha_root) - lehat(i,1)*lattice.N(i,1)*cos(state.betha) + lehat(i,3)*lattice.N(i,3)*cos(state.betha) + lehat(i,1)*lattice.N(i,2)*sin(state.betha) + 2*lehat(i,1)*lattice.N(i,1)*cos(state.alpha_root)^2*cos(state.betha) - 2*lehat(i,3)*lattice.N(i,3)*cos(state.alpha_root)^2*cos(state.betha) - 2*lehat(i,1)*lattice.N(i,2)*cos(state.alpha_root)^2*sin(state.betha) + 2*lehat(i,3)*lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha)^2*sin(state.alpha_root) - 2*lehat(i,3)*lattice.N(i,2)*cos(state.alpha_root)*cos(state.betha)*sin(state.alpha_root)*sin(state.betha),...
%         -((lehat(i,2)*cos(state.betha) + lehat(i,1)*sin(state.betha))*(lattice.N(i,3) - 2*lattice.N(i,3)*cos(state.alpha_root)^2 + 2*lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha)*sin(state.alpha_root) - 2*lattice.N(i,2)*cos(state.alpha_root)*sin(state.alpha_root)*sin(state.betha)))];
results.dF_ar(:,:) = state.AS(i)^2*state.rho*[invW2(:,i).*(-(lehat(i,3)*(cos(state.alpha_root)*sin(state.betha) + DWY*lehat(i,2)*invW2(:,i).*...
    (lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha))...
    - lehat(i,3)*(sin(state.alpha_root)*sin(state.betha) - DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3) *cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root)...
    + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha))...
    + lehat(i,2)*(cos(state.alpha_root) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))...
    *(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha))...
    + lehat(i,2)*(sin(state.alpha_root) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))),...
    ...
    invW2(:,i).*((lehat(i,3)*(cos(state.betha)*sin(state.alpha_root) + DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))...
    *(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)) - lehat(i,3)*(cos(state.alpha_root)*cos(state.betha)...
    - DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha_root)...
    - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)) + lehat(i,1)*(cos(state.alpha_root) + DWZ*lehat(i,3)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha_root)...
    - lehat(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha)...
    - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)) + lattice.N(i,1)*(sin(state.alpha_root) + DWZ*lattice.N(i,3)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha)...
    - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))),...
    ...
    invW2(:,i).*((lehat(i,2)*(cos(state.alpha_root)*cos(state.betha) - DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))...
    *(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)) - lehat(i,2)*(cos(state.betha)*sin(state.alpha_root)...
    + DWX*lehat(i,1)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root) + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha_root)...
    + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)) + lehat(i,1)*(cos(state.alpha_root)*sin(state.betha) + DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3)*sin(state.alpha_root)...
    + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root)...
    + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)) - lattice.N(i,1)*(sin(state.alpha_root)*sin(state.betha) - DWY*lehat(i,2)*invW2(:,i).*(lattice.N(i,3)*cos(state.alpha_root) - lattice.N(i,1)*cos(state.betha)*sin(state.alpha_root)...
    + lattice.N(i,2)*sin(state.alpha_root)*sin(state.betha)))*(lattice.N(i,3)*sin(state.alpha_root) + lattice.N(i,1)*cos(state.alpha_root)*cos(state.betha) - lattice.N(i,2)*cos(state.alpha_root)*sin(state.betha))))];

results.dM_ar = cross(squeeze(C3(:,1,:)),results.dF_ar);
results.dFORCE_ar = sum(results.dF_ar,1);
results.NP = [0,results.dFORCE_ar(3),-results.dFORCE_ar(2);...
    -results.dFORCE_ar(3),0,results.dFORCE_ar(1);...
    results.dFORCE_ar(2),results.dFORCE_ar(1),0]\sum(cross(lattice.COLLOC,results.dF_ar))';
results.dMOMENTS_ar = sum(results.dM_ar,1);

results.F=F;
results.FORCE=sum(F,1);						%Total force
M=cross(C3,F,3);			                 %Moments per panel
results.M=M;
results.MOMENTS=sum(M,1);					%Summing up moments	
results.gamma=gamma;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NEW DOWNWASH FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dw,DW]=fastdw(lattice)
one_by_four_pi=1/(4*pi);

[psize vsize void]=size(lattice.VORTEX);


%disp('running right')
%psize=size(lattice.COLLOC,1);
lemma=ones(1,psize);

LDW=zeros(psize,psize,7,3);

mCOLLOC(:,:,1)=lattice.COLLOC(:,1)*lemma;
mCOLLOC(:,:,2)=lattice.COLLOC(:,2)*lemma;
mCOLLOC(:,:,3)=lattice.COLLOC(:,3)*lemma;

mN(:,:,1)=lattice.N(:,1)*lemma;
mN(:,:,2)=lattice.N(:,2)*lemma;
mN(:,:,3)=lattice.N(:,3)*lemma;

for j=1:(vsize-1)
    
    lr1(:,:,1)=(lattice.VORTEX(:,j,1)*lemma)';
    lr1(:,:,2)=(lattice.VORTEX(:,j,2)*lemma)';
    lr1(:,:,3)=(lattice.VORTEX(:,j,3)*lemma)';
    
    lr2(:,:,1)=(lattice.VORTEX(:,j+1,1)*lemma)';
    lr2(:,:,2)=(lattice.VORTEX(:,j+1,2)*lemma)';
    lr2(:,:,3)=(lattice.VORTEX(:,j+1,3)*lemma)'; 
    
    r1=lr1-mCOLLOC;
    r2=lr2-mCOLLOC;
    
    warning off
    LDW(:,:,j,:)=mega(r1,r2);
    warning on
end
LDW(find((isnan(LDW(:,:,:,:)))))=0;

DW=-squeeze(sum(LDW,3))*one_by_four_pi;

dw=sum(DW.*mN,3);
end


function[DW2]=mega(r1,r2)
%% First part
F1=cross(r1,r2,3);

LF1=(sum(F1.^2,3));



F2(:,:,1)=F1(:,:,1)./(LF1);
F2(:,:,2)=F1(:,:,2)./(LF1);
F2(:,:,3)=F1(:,:,3)./(LF1);
%clear('F1')


%% Next part

Lr1=sqrt(sum(r1.^2,3)); 
Lr2=sqrt(sum(r2.^2,3));


R1(:,:,1)=r1(:,:,1)./Lr1;
R1(:,:,2)=r1(:,:,2)./Lr1;
R1(:,:,3)=r1(:,:,3)./Lr1;

R2(:,:,1)=r2(:,:,1)./Lr2;
R2(:,:,2)=r2(:,:,2)./Lr2;
R2(:,:,3)=r2(:,:,3)./Lr2;



L1=(R2-R1);
%clear('R1','R2')



%% Third part
R0=(r2-r1);

radial_distance=sqrt((LF1./(sum(R0.^2,3))));


%% combinging 2 and 3
L2=  R0(:,:,1).*L1(:,:,1)...
    +R0(:,:,2).*L1(:,:,2)...
    +R0(:,:,3).*L1(:,:,3);



%% Downwash
DW(:,:,1)=F2(:,:,1).*L2;
DW(:,:,2)=F2(:,:,2).*L2;
DW(:,:,3)=F2(:,:,3).*L2;

near=config('near');

DW2(:,:,1)=DW(:,:,1).*(1-(radial_distance<near));
DW2(:,:,2)=DW(:,:,2).*(1-(radial_distance<near));
DW2(:,:,3)=DW(:,:,3).*(1-(radial_distance<near));


end