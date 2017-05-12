%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lattice] = generateLattice(geo,S_ref,C_ref,B_ref,mac_pos,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generateLattice: Function for Dynamic TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the lattice structure for the calculation of
% the vortex lattice method. It convernts the finite element states
% into a usable lattice.
% NOTE: Current calculations are not for wings with sweep,
%       taper, or other variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Function for Dynamic TORNADO Aeroelastic Simulator					
% Called by:	---------FILL IN LATER----------												
% Calls:	normals4
%           wakesetup2
%           MATLAB 5.2 std fcns
% Inputs:   geo - Geometry of the plane
%                   ref_point - aircraft geometry refrence point
%                   CG - aircraft center of gravity
%                   nwing - number of wings
%                   symetric - is wing symetric about y axis
%                   Wings - Array of wing structures
%                           wing - Wing structure
%                                   meanCord - mean cord line
%                                   cord - cord length
%                                   span - wing span
%                                   L - Length array
%                                   SegNum - number of spanwise segments
%                                   cordNum - number of cordwise segments
%                                   Theta - wing twist array
%                                   Z - z displacement array
%                                   PhiZ - Z slope array
%                                   X - x displacement array
%                                   PhiX - x slope array
%                                   start - wing start location
%                                   Controls - array of control surfaces
%                                               surface - corners of the
%                                               control surface area
%                                               hinge - hinge location,
%                                               [x1,y1,x2,y2]
%                                               rotation - rotation degree
%                                               following right hand rule
%           numCordB - Number of lattice partitions in the cord of the body
%           bw - Body Width
%           bc - Body Length
%           numSpanB - Number of partitions in the span direction of the
%           body
%           state - Structure containing simulation state:
%                       alpha - angle of attack
%                       betha - angle of side sweep
%                       P - Roll angular velocity [deg/s], set to zero
%                       Q - Pitch angular velocity [deg/s], set to zero
%                       R - Yaw angular velocity [deg/s], set to zero
%                       adot - Angle of attack time derivative, (Alpha_dot),
%                           [deg/s], set to zero
%                       bdot - Angle of sideslip time derivative, (Beta_dot),
%                           [deg/s], set to zero
%                       AS - Airspeed
%                       rho - Air density
%                       ALT - altitude, set to zero because air density is
%                          provided
%                       alphadot - likely legacy code, set to zero
%                       bethadot - likely legacy code, set to zero
%                       pgcorr - Apply Prandtl-Glauert Correction [0 1], zero
%                          is false, one is true, set to zero
% Output:   lattice - lattice structure containing:
%                   COLLOC - cordNum*(SegNum+numSpanB)by 3 positional
%                   matrix of collocation points
%                   VORTEX - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                   containing x,y,z positions for each corner point of the
%                   vortex
%                   N - Normal vectors of VORTEX
%                   XYZ - cordNum*(SegNum+numSpanB)by 4 of 3 matrices
%                   containing x,y,z positions for each external corner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intializing lattice
lattice.COLLOC=[];
lattice.VORTEX=[];
lattice.N=[];
lattice.XYZ=[];

%Initalize Position arrays
x1n = [];
x2n = [];
xTrailR = [];
xTrailL = [];
xLeadR = [];
xLeadL = [];
zTrailR = [];
zLeadR = [];
zTrailL = [];
zLeadL = [];
yTrailR = [];
yLeadR = [];
yTrailL = [];
yLeadL = [];
y1n = [];
y2n = [];
z1n = [];
z2n = [];
ym = [];
xm = [];
zm = [];

%Iterate through all the wings 
for k = 1:geo.nwing
    for j = 1:geo.Wings(k).wing.SegNum
        %Interpolate the the mean camber for the number of cordwise
        %sections
        camberX = [0:geo.Wings(k).wing.chord/geo.Wings(k).wing.cordNum:geo.Wings(k).wing.chord];
        camberY = interp1(geo.Wings(k).wing.meanCamber(j).camber(:,1)*geo.Wings(k).wing.chord(j),geo.Wings(k).wing.meanCamber(j).camber(:,2)*geo.Wings(k).wing.chord(j),camberX,'pchip');
        %Go through each cord section
        for i = 1:geo.Wings(k).wing.cordNum
            %Create the chord of segment wing slice on the mean camber line
            SegCord = [camberX(i+1)-camberX(i),camberY(i+1)-camberY(i)];
            %Calculate quarter cord, 3/4 cord and leading and trealing edge
            %of the panel segment
            quartCord = .25*SegCord+[camberX(i),camberY(i)];
            threeQuartCord = .75*SegCord+[camberX(i),camberY(i)];
            trail = [camberX(i+1),camberY(i+1)];
            lead = [camberX(i),camberY(i)];
            
            %Calculate the collocation point
            alpha = 0.5*geo.Wings(k).wing.Theta(j+1)+0.5*geo.Wings(k).wing.Theta(j);
            collocationCord = [interp1(geo.Wings(k).wing.Y,SegCord(1),0.5*geo.Wings(k).wing.Y(j+1)+0.5*geo.Wings(k).wing.Y(j));...
                interp1(geo.Wings(k).wing.Y,SegCord(2),0.5*geo.Wings(k).wing.Y(j+1)+0.5*geo.Wings(k).wing.Y(j))];
            if j < (wing.SegNum-1)/2+1
                zm = [zm;0.5*geo.Wings(k).wing.Z(j+1)+.125*geo.Wings(k).wing.PhiZ(j+1)+0.5*geo.Wings(k).wing.Z(j)-0.125*geo.Wings(k).wing.PhiZ(j)+[sin(alpha),cos(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+geo.Wings(k).wing.start(3)];
                xm = [xm;0.5*geo.Wings(k).wing.X(j+1)+.125*geo.Wings(k).wing.PhiX(j+1)+0.5*geo.Wings(k).wing.X(j)-0.125*geo.Wings(k).wing.PhiX(j)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+geo.Wings(k).wing.start(1)];
            else
                zm = [zm;0.5*geo.Wings(k).wing.Z(j)+.125*geo.Wings(k).wing.PhiZ(j)+0.5*geo.Wings(k).wing.Z(j+1)-0.125*geo.Wings(k).wing.PhiZ(j+1)+[sin(alpha),cos(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+geo.Wings(k).wing.start(3)];
                xm = [xm;0.5*geo.Wings(k).wing.X(j)+.125*geo.Wings(k).wing.PhiX(j)+0.5*geo.Wings(k).wing.X(j+1)-0.125*geo.Wings(k).wing.PhiX(j+1)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+geo.Wings(k).wing.start(1)];
            end
            ym = [ym;(geo.Wings(k).wing.Y(j)+geo.Wings(k).wing.Y(j+1))/2];
            %Calculate X and Z corrdinated based off twist
            x1n = [x1n;geo.Wings(k).wing.X(j)+diag([cos(geo.Wings(k).wing.Theta(j)),-sin(geo.Wings(k).wing.Theta(j))]*quartCord(:,j))+geo.Wings(k).wing.start(1)];
            x2n = [x2n;geo.Wings(k).wing.X(j+1)+diag([cos(geo.Wings(k).wing.Theta(j+1)),-sin(geo.Wings(k).wing.Theta(j+1))]*quartCord(:,j+1))+geo.Wings(k).wing.start(1)];
            xTrailR = [xTrailR;geo.Wings(k).wing.X(j)+diag([cos(geo.Wings(k).wing.Theta(j)),-sin(geo.Wings(k).wing.Theta(j))]*trail(:,j))+geo.Wings(k).wing.start(1)];
            xTrailL = [xTrailL;geo.Wings(k).wing.X(j+1)+diag([cos(geo.Wings(k).wing.Theta(j+1)),-sin(geo.Wings(k).wing.Theta(j+1))]*trail(:,j+1))+geo.Wings(k).wing.start(1)];
            xLeadR = [xLeadR;geo.Wings(k).wing.X(j)+diag([cos(geo.Wings(k).wing.Theta(j)),-sin(geo.Wings(k).wing.Theta(j))]*lead(:,j))+geo.Wings(k).wing.start(1)];
            xLeadL = [xLeadL;geo.Wings(k).wing.X(j+1)+diag([cos(geo.Wings(k).wing.Theta(j+1)),-sin(geo.Wings(k).wing.Theta(j+1))]*lead(:,j+1))+geo.Wings(k).wing.start(1)];
            %
            z1n = [z1n;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j)),cos(geo.Wings(k).wing.Theta(j))]*quartCord(:,j))+geo.Wings(k).wing.start(3)];
            z2n = [z2n;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j+1)),cos(geo.Wings(k).wing.Theta(j+1))]*quartCord(:,j+1))+geo.Wings(k).wing.start(3)];
            zTrailR = [zTrailR;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j)),cos(geo.Wings(k).wing.Theta(j))]*trail(:,j))+geo.Wings(k).wing.start(3)];
            zTrailL = [zTrailL;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j+1)),cos(geo.Wings(k).wing.Theta(j+1))]*trail(:,j+1))+geo.Wings(k).wing.start(3)];
            zLeadR = [zLeadR;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j)),cos(geo.Wings(k).wing.Theta(j))]*lead(:,j))+geo.Wings(k).wing.start(3)];
            zLeadL = [zLeadL;geo.Wings(k).wing.Z(j)+diag([sin(geo.Wings(k).wing.Theta(j+1)),cos(geo.Wings(k).wing.Theta(j+1))]*lead(:,j+1))+geo.Wings(k).wing.start(3)];
            %
            yTrailR = [yTrailR;geo.Wings(k).wing.Y(1:end-1)+geo.Wings(k).wing.start(2)];
            yTrailL = [yTrailL;geo.Wings(k).wing.Y(2:end)+geo.Wings(k).wing.start(2)];
            yLeadR = [yLeadR;geo.Wings(k).wing.Y(1:end-1)+geo.Wings(k).wing.start(2)];
            yLeadL = [yLeadL;geo.Wings(k).wing.Y(2:end)+geo.Wings(k).wing.start(2)];
            y1n = [y1n;geo.Wings(k).wing.Y(1:end-1)+geo.Wings(k).wing.start(2)];
            y2n = [y2n;geo.Wings(k).wing.Y(2:end)+geo.Wings(k).wing.start(2)];
        end
    end
end

lattice.COLLOC = [xm,ym,zm];

lattice.VORTEX(:,:,1) = [xTrailR,x1n,x2n,xTrailR];
lattice.VORTEX(:,:,2) = [y1n,y1n,y2n,y2n];
lattice.VORTEX(:,:,3) = [z1n,z1n,z2n,z2n];

lattice.XYZ(:,:,1) = [xLeadR,xLeadL,xTrailL,xTrailR,xLeadR];
lattice.XYZ(:,:,2) = [yLeadR,yLeadL,yTrailL,yTrailR,yTrailR];
lattice.XYZ(:,:,3) = [zLeadR,zLeadL,zTrailL,zTrailR,zLeadR];

S = zeros(1,geo.nx*geo.ny);
lattice.N=normals4(lattice.COLLOC,lattice.VORTEX,S);

[ref]=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos);
if state.AS~=0   %appending wake lattice points (farpoints)
    lattice=wakesetup2(lattice,state,ref); %setting up wake legs.
    temporary=lattice.VORTEX(:,[1 3 4 6],:);
    temporary(:,1,3)=temporary(:,2,3); %Flattening wake
    temporary(:,4,3)=temporary(:,3,3); %Flattening wake
    lattice.VORTEX=temporary;
else
    terror(13)
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [normal]=normals4(colloc,vortex,C_Slope)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALS: Essential function for TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the normals to								
% each panel. Two vectors in the plane, the ones between	
% the vortex points and the collocation point, defines	
% the panel plane. Together with the vortex orientation	
% the orientation of the normal is defined.					
% Output normals are normalized.									
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautics	
%				copyright 2000											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Essential function for TORNADO					
% Called by:	setup												
% Calls:			trot												
%					MATLAB 5.2 std fcns							
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=[];
step=size(colloc);
[d e f]=size(vortex);
a=e/2;
b=a+1;
      
for t=1:step	%Looping through panels
   alpha=C_Slope(t);
   
   	for s=1:3						%Looping Through Dimensions.
      	ra(s)=vortex(t,a,s);
      	rb(s)=vortex(t,b,s);
      	rc(s)=colloc(t,s);
      end
        r0=rb-ra;
        r0(1)=0;                    %fix to get normals to not point the right way
      	r1=rc-ra;
      	r2=rc-rb;
   		n=cross(r1,r2);				%Passus to determine normal
      	nl=sqrt(sum((n.^2),2));    %of panel at collocationpoint.
    		R=n/nl;							%Normalizing normal.
         R2=trot3(r0,R,-alpha);		%rotating wha trot
         N=[N;R2']; 
end

normal=N;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[p2]=trot3(hinge,p,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TROT: Auxillary rotation function			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotates point p around hinge alpha rads.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ref: 	Råde, Westergren, BETA 4th ed,   
%			studentlitteratur, 1998			    	
%			pp:107-108							   	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: 	Tomas Melin, KTH,Department of%
% 				aeronautics, Copyright 2000	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Context:	Auxillary function for			
%				TORNADO.								
% Called by: setrudder, normals			
% Calls:		norm (MATLAB std fcn)			
%				sin			"						
%				cos			"						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP:		Hinge=vector around rotation  
%						takes place.				
%				p=point to be rotated			
%				alpha=radians of rotation		
%				3D-workspace						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=hinge(1);
b=hinge(2);
c=hinge(3);

rho=sqrt(a^2+b^2);
r=sqrt(a^2+b^2+c^2);

if r==0
   cost=0;
   sint=1;
else
   cost=c/r;
   sint=rho/r;
end

if rho==0
   cosf=0;
   sinf=1;
else
   cosf=a/rho;
	sinf=b/rho;
end   

cosa=cos(alpha);
sina=sin(alpha);

RZF=[[cosf -sinf 0];[sinf cosf 0];[0 0 1]];
RYT=[[cost 0 sint];[0 1 0];[-sint 0 cost]];
RZA=[[cosa -sina 0];[sina cosa 0];[0 0 1]];
RYMT=[[cost 0 -sint];[0 1 0];[sint 0 cost]];
RZMF=[[cosf sinf 0];[-sinf cosf 0];[0 0 1]];

P=RZF*RYT*RZA*RYMT*RZMF;
p2=P*p';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lattice]=wakesetup2(lattice,state,ref)
infdist=config('infinity');
if isempty(infdist)
	infdist=6*ref.b_ref;
end

[a,b,c]=size(lattice.VORTEX);
V2=lattice.VORTEX;
c=[1 b];

infx=infdist*cos(state.alpha_root)*cos(state.betha);
infy=-infdist*sin(state.betha);
infz=infdist*sin(state.alpha_root)*cos(state.betha);

for t=1:a
   for s=1:2
   	  x=infx+lattice.VORTEX(t,c(s),1);
      y=infy+lattice.VORTEX(t,c(s),2);
      z=infz+lattice.VORTEX(t,c(s),3);
      
      psi=state.P/state.U_inf*x;
  	  theta=state.Q/state.U_inf*x;
   	  fi=state.R/state.U_inf*x;
      
      dx(t,s)=-x*(2-cos(theta)-cos(fi));
   	  dy(t,s)=+sin(psi)*z-sin(fi)*x+(1-cos(psi))*y;
      dz(t,s)=sin(theta)*x-sin(psi)*y+(1-cos(psi))*z;
      
   end
end

for i=1:a
   INF1(i,1,1)=lattice.VORTEX(i,1,1)+infx+dx(i,1);
   INF1(i,1,2)=lattice.VORTEX(i,1,2)+infy+dy(i,1);
   INF1(i,1,3)=lattice.VORTEX(i,1,3)+infz+dz(i,1);
   
   INF2(i,1,1)=lattice.VORTEX(i,b,1)+infx+dx(i,2);
   INF2(i,1,2)=lattice.VORTEX(i,b,2)+infy+dy(i,2);
   INF2(i,1,3)=lattice.VORTEX(i,b,3)+infz+dz(i,2);
end

lattice.VORTEX=[INF1(:,1,:) V2(:,:,:) INF2(:,1,:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xm,ym,zm,x1n,x2n,xTrailR,xTrailL,xLeadR,xLeadL,y1n,y2n,yTrailR,yTrailL,yLeadR,yLeadL,z1n,z2n,zTrailR,zTrailL,zLeadR,zLeadL] = ...
    spanWisePopulate(wing,SegCord,quartCord,threeQuartCord,trail,lead,xm,ym,zm,x1n,x2n,xTrailR,xTrailL,xLeadR,xLeadL,y1n,y2n,yTrailR,yTrailL,yLeadR,yLeadL,z1n,z2n,zTrailR,zTrailL,zLeadR,zLeadL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spanWisePopulate: polulates components for
%   spanwise lattice section of wing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: 	Nick Cramer, UCSC,Department of%
% 			Computer Engineering, Copyright 2015	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Context:	Auxillary function for			
%				Dynamic TORNADO.								
% Called by: generatelattice		
% Calls:		diag (MATLAB std fcn)								
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:    wing - Wing structure
%                                   meanCord - mean cord line
%                                   cord - cord length
%                                   span - wing span
%                                   L - Length array
%                                   SegNum - number of spanwise segments
%                                   cordNum - number of cordwise segments
%                                   Theta - wing twist array
%                                   Z - z displacement array
%                                   PhiZ - Z slope array
%                                   X - x displacement array
%                                   PhiX - x slope array
%                                   start - wing start location
%                                   Controls - array of control surfaces
%                                               surface - corners of the
%                                               control surface area
%                                               hinge - hinge location,
%                                               [x1,y1,x2,y2]
%                                               rotation - rotation degree
%                                               following right hand rule
%           quartCord - array of quarter cord coordninates
%           threeQuartCord - array of three quarter cord coordinates
%           trail - array of the trailing edge coordinates
%           lead - array of the leading edge coordinates
%           xm - array of the x component for the collication point
%           ym - array of the y component for the collication point
%           zm - array of the z component for the collication point
%           x1n - array of the x component for the left most point of the
%               horseshoe
%           x2n - array of the x component for the right most point of the
%               horseshoe
%           xTrailR - array of the x component for the right most trailing
%               edge of the panel
%           xTrailL - array of the x component for the left most trailing
%               edge of the panel
%           xLeadR - array of the x component for the right most leading
%               edge of the panel
%           xLeadL - array of the x component for the left most trailing
%               edge of the panel
%           y1n - array of the y component for the left most point of the
%               horseshoe
%           y2n - array of the y component for the right most point of the
%               horseshoe
%           yTrailR - array of the y component for the right most trailing
%               edge of the panel
%           yTrailL - array of the y component for the right most trailing
%               edge of the panel
%           yLeadR - array of the y component for the right most leading
%               edge of the panel
%           yLeadL - array of the y component for the left most trailing
%               edge of the panel
%           z1n - array of the z component for the left most point of the
%               horseshoe
%           z2n - array of the z component for the right most point of the
%               horseshoe
%           zTrailR - array of the z component for the right most trailing
%               edge of the panel
%           zTrailL - array of the z component for the right most trailing
%               edge of the panel
%           zLeadR - array of the z component for the right most leading
%               edge of the panel
%           zLeadL - array of the z component for the left most trailing
%               edge of the panel
% Outputs:  xm - array of the x component for the collication point
%           ym - array of the y component for the collication point
%           zm - array of the z component for the collication point
%           x1n - array of the x component for the left most point of the
%               horseshoe
%           x2n - array of the x component for the right most point of the
%               horseshoe
%           xTrailR - array of the x component for the right most trailing
%               edge of the panel
%           xTrailL - array of the x component for the left most trailing
%               edge of the panel
%           xLeadR - array of the x component for the right most leading
%               edge of the panel
%           xLeadL - array of the x component for the left most trailing
%               edge of the panel
%           y1n - array of the y component for the left most point of the
%               horseshoe
%           y2n - array of the y component for the right most point of the
%               horseshoe
%           yTrailR - array of the y component for the right most trailing
%               edge of the panel
%           yTrailL - array of the y component for the right most trailing
%               edge of the panel
%           yLeadR - array of the y component for the right most leading
%               edge of the panel
%           yLeadL - array of the y component for the left most trailing
%               edge of the panel
%           z1n - array of the z component for the left most point of the
%               horseshoe
%           z2n - array of the z component for the right most point of the
%               horseshoe
%           zTrailR - array of the z component for the right most trailing
%               edge of the panel
%           zTrailL - array of the z component for the right most trailing
%               edge of the panel
%           zLeadR - array of the z component for the right most leading
%               edge of the panel
%           zLeadL - array of the z component for the left most trailing
%               edge of the panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Update control sufrace part of the code to be the same calculation of
%colocation points as the non controls part!!!

%Check to see if there are any control surfaces
if isempty(wing.Controls)
    switch wing.plane
        case 1
            for j = 1:wing.SegNum
                %Calculate the collocation point
                alpha = 0.5*wing.Theta(j+1)+0.5*wing.Theta(j);
                collocationCord = [interp1(wing.Y,SegCord(1,:),0.5*wing.Y(j+1)+0.5*wing.Y(j));...
                    interp1(wing.Y,SegCord(2,:),0.5*wing.Y(j+1)+0.5*wing.Y(j))];
                if j < (wing.SegNum-1)/2+1
                    zm = [zm;0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j)+[sin(alpha),cos(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(3)];
                    xm = [xm;0.5*wing.X(j+1)+.125*wing.PhiX(j+1)+0.5*wing.X(j)-0.125*wing.PhiX(j)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(1)];
                else
                    zm = [zm;0.5*wing.Z(j)+.125*wing.PhiZ(j)+0.5*wing.Z(j+1)-0.125*wing.PhiZ(j+1)+[sin(alpha),cos(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(3)];
                    xm = [xm;0.5*wing.X(j)+.125*wing.PhiX(j)+0.5*wing.X(j+1)-0.125*wing.PhiX(j+1)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(1)];
                end
                ym = [ym;(wing.Y(j)+wing.Y(j+1))/2];
                %Calculate X and Z corrdinated based off twist
                x1n = [x1n;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*quartCord(:,j))+wing.start(1)];
                x2n = [x2n;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(1)];
                xTrailR = [xTrailR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*trail(:,j))+wing.start(1)];
                xTrailL = [xTrailL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*trail(:,j+1))+wing.start(1)];
                xLeadR = [xLeadR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*lead(:,j))+wing.start(1)];
                xLeadL = [xLeadL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*lead(:,j+1))+wing.start(1)];
                %
                z1n = [z1n;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*quartCord(:,j))+wing.start(3)];
                z2n = [z2n;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(3)];
                zTrailR = [zTrailR;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*trail(:,j))+wing.start(3)];
                zTrailL = [zTrailL;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*trail(:,j+1))+wing.start(3)];
                zLeadR = [zLeadR;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*lead(:,j))+wing.start(3)];
                zLeadL = [zLeadL;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*lead(:,j+1))+wing.start(3)];
            end
            %
            yTrailR = [yTrailR;wing.Y(1:end-1)+wing.start(2)];
            yTrailL = [yTrailL;wing.Y(2:end)+wing.start(2)];
            yLeadR = [yLeadR;wing.Y(1:end-1)+wing.start(2)];
            yLeadL = [yLeadL;wing.Y(2:end)+wing.start(2)];
            y1n = [y1n;wing.Y(1:end-1)+wing.start(2)];
            y2n = [y2n;wing.Y(2:end)+wing.start(2)];
        case 2
            for j = 1:wing.SegNum
                %Calculate the collocation point
                alpha = 0.5*wing.Theta(j+1)+0.5*wing.Theta(j);
                collocationCord = [interp1(wing.Z,SegCord(1,:),0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j));...
                    interp1(wing.Z,SegCord(2,:),0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j))];
                if j < (wing.SegNum-1)/2+1
                    zm = [zm;0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j)+wing.start(3)];
                    xm = [xm;0.5*wing.X(j+1)+.125*wing.PhiX(j+1)+0.5*wing.X(j)-0.125*wing.PhiX(j)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(1)];
                else
                    zm = [zm;0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j)+wing.start(3)];
                    xm = [xm;0.5*wing.X(j)+.125*wing.PhiX(j)+0.5*wing.X(j+1)-0.125*wing.PhiX(j+1)+[cos(alpha),-sin(alpha)]*(.75*collocationCord+(lead(:,j)+lead(:,j+1))/2)+wing.start(1)];
                end
                ym = [ym;(wing.Y(j)+wing.Y(j+1))/2+[sin(alpha),cos(alpha)]*threeQuartCord(:,j)+wing.start(2)];
                %Calculate X and Z corrdinated based off twist
                x1n = [x1n;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*quartCord(:,j))+wing.start(1)];
                x2n = [x2n;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(1)];
                xTrailR = [xTrailR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*trail(:,j))+wing.start(1)];
                xTrailL = [xTrailL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*trail(:,j+1))+wing.start(1)];
                xLeadR = [xLeadR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*lead(:,j))+wing.start(1)];
                xLeadL = [xLeadL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*lead(:,j+1))+wing.start(1)];
                %
                y1n = [y1n;wing.Y(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*quartCord(:,j))+wing.start(2)];
                y2n = [y2n;wing.Y(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(2)];
                yTrailR = [yTrailR;wing.Y(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*trail(:,j))+wing.start(2)];
                yTrailL = [yTrailL;wing.Y(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*trail(:,j+1))+wing.start(2)];
                yLeadR = [yLeadR;wing.Y(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*lead(:,j))+wing.start(2)];
                yLeadL = [yLeadL;wing.Y(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*lead(:,j+1))+wing.start(2)];
            end
            %
            z1n = [z1n;wing.Z(1:end-1)+wing.start(3)];
            z2n = [z2n;wing.Z(2:end)+wing.start(3)];
            zTrailR = [zTrailR;wing.Z(1:end-1)+wing.start(3)];
            zTrailL = [zTrailL;wing.Z(2:end)+wing.start(3)];
            zLeadR = [zLeadR;wing.Z(1:end-1)+wing.start(3)];
            zLeadL = [zLeadL;wing.Z(2:end)+wing.start(3)];
    end
else
    %Calculate the y coordinates for panels
    y1nTemp = wing.Y(1:end-1);
    y2nTemp = wing.Y(2:end);

    for j = 1:wing.SegNum
        %Calculate x position leading vortex edges
        x1nTemp = wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*quartCord(:,j))+wing.start(1);
        x2nTemp = wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(1);
        
        %Determine if current pannel is within a control surface
        isControl = zeros(1,length(wing.Controls));
        for k = 1:length(wing.Controls)
            isControl(k) = ((x1nTemp>wing.Controls(k).surface(1,1)&x1nTemp<wing.Controls(k).surface(2,1))|(x2nTemp>wing.Controls(k).surface(1,1)&x2nTemp<wing.Controls(k).surface(2,1)))&...
                ((y1nTemp(j)>wing.Controls(k).surface(1,2)&y1nTemp(j)<wing.Controls(k).surface(2,2))|(y2nTemp(j)>wing.Controls(k).surface(1,2)&y2nTemp(j)<wing.Controls(k).surface(2,2)));
        end
        
        %Calculate the collocation point
        alpha = 0.5*wing.Theta(j+1)+0.5*wing.Theta(j);
        if j < (wing.SegNum-1)/2+1
            zm = [zm;0.5*wing.Z(j+1)+.125*wing.PhiZ(j+1)+0.5*wing.Z(j)-0.125*wing.PhiZ(j)+[sin(alpha),cos(alpha)]*threeQuartCord(:,j)+wing.start(3)];
            xm = [xm;0.5*wing.X(j+1)+.125*wing.PhiX(j+1)+0.5*wing.X(j)-0.125*wing.PhiX(j)+[cos(alpha),-sin(alpha)]*threeQuartCord(:,j)+wing.start(1)];
        else
            zm = [zm;0.5*wing.Z(j)+.125*wing.PhiZ(j)+0.5*wing.Z(j+1)-0.125*wing.PhiZ(j+1)+[sin(alpha),cos(alpha)]*threeQuartCord(:,j)+wing.start(3)];
            xm = [xm;0.5*wing.X(j)+.125*wing.PhiX(j)+0.5*wing.X(j+1)-0.125*wing.PhiX(j+1)+[cos(alpha),-sin(alpha)]*threeQuartCord(:,j)+wing.start(1)];
        end
        ym = [ym;(wing.Y(j)+wing.Y(j+1))/2];
        
        %Calculate X and Z corrdinated based off twist
        x1n = [x1n;x1nTemp];
        x2n = [x2n;x2nTemp];
        xTrailR = [xTrailR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*trail(:,j))+wing.start(1)];
        xTrailL = [xTrailL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*trail(:,j+1))+wing.start(1)];
        xLeadR = [xLeadR;wing.X(j)+diag([cos(wing.Theta(j)),-sin(wing.Theta(j))]*lead(:,j))+wing.start(1)];
        xLeadL = [xLeadL;wing.X(j+1)+diag([cos(wing.Theta(j+1)),-sin(wing.Theta(j+1))]*lead(:,j+1))+wing.start(1)];
        %
        z1n = [z1n;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*quartCord(:,j))+wing.start(3)];
        z2n = [z2n;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*quartCord(:,j+1))+wing.start(3)];
        zTrailR = [zTrailR;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*trail(:,j))+wing.start(3)];
        zTrailL = [zTrailL;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*trail(:,j+1))+wing.start(3)];
        zLeadR = [zLeadR;wing.Z(j)+diag([sin(wing.Theta(j)),cos(wing.Theta(j))]*lead(:,j))+wing.start(3)];
        zLeadL = [zLeadL;wing.Z(j)+diag([sin(wing.Theta(j+1)),cos(wing.Theta(j+1))]*lead(:,j+1))+wing.start(3)];
        if any(isControl)
            %Determine which control surface the point is contained in
            cSurface = find(isControl);
            
            %Get angle that the control surface is rotated
            angle = wing.Controls(cSurface).rotation*pi/180;
            
            %Generate geometric constants
            a = wing.Controls(cSurface).hinge(1,1);
            b = wing.Controls(cSurface).hinge(1,2);
            c = wing.Controls(cSurface).hinge(1,3);
            u = wing.Controls(cSurface).hinge(2,1)-wing.Controls(cSurface).hinge(1,1);
            v = wing.Controls(cSurface).hinge(2,2)-wing.Controls(cSurface).hinge(1,2);
            w = wing.Controls(cSurface).hinge(2,3)-wing.Controls(cSurface).hinge(1,3);
            L = norm([u,v,w]);
            alpha = atan2(u,v);
            beta = atan2(w,v);
            
            %Create Rotation Matrix about the hinge vector
            Tabc = [eye(3),[-a;-b;-c];[0,0,0,1]];
            Rz = [cos(alpha),-sin(alpha),0,0;sin(alpha),cos(alpha),0,0;zeros(2),eye(2)];
            Rx = [1,0,0,0;0,cos(beta),-sin(beta),0;0,sin(beta),cos(beta),0;0,0,0,1];
            Ry = [cos(angle),0,sin(angle),0;0,1,0,0;-sin(angle),0,cos(angle),0;0,0,0,1];
            iRx = [1,0,0,0;0,cos(-beta),-sin(-beta),0;0,sin(-beta),cos(-beta),0;0,0,0,1];
            iRz = [cos(-alpha),-sin(-alpha),0,0;sin(-alpha),cos(-alpha),0,0;zeros(2),eye(2)];
            iTabc = [eye(3),[a;b;c];[0,0,0,1]];
            rotationMatrix = iTabc*iRz*iRx*Ry*Rx*Rz*Tabc;
            
            %Rotate quarter cord and edge points
            temp1n = rotationMatrix*[x1n(end);y1nTemp(j);z1n(end);1];
            temp2n = rotationMatrix*[x2n(end);y2nTemp(j);z2n(end);1];
            tempTrailR = rotationMatrix*[xTrailR(end);y1nTemp(j);zTrailR(end);1];
            tempTrailL = rotationMatrix*[xTrailL(end);y2nTemp(j);zTrailL(end);1];
            tempLeadR = rotationMatrix*[xLeadR(end);y1nTemp(j);zLeadR(end);1];
            tempLeadL = rotationMatrix*[xLeadL(end);y2nTemp(j);zLeadL(end);1];
            tempM = rotationMatrix*[xm(end);ym(end);zm(end);1];
            
            %Save rotated points into position array
            x1n(end) = temp1n(1);
            y1nTemp(j) = temp1n(2);
            z1n(end) = temp1n(3);
            x2n(end) = temp2n(1);
            y2nTemp(j) = temp2n(2);
            z2n(end) = temp2n(3);
            xTrailR(end) = tempTrailR(1);
            %Note: Not an issue right now but there will need to evenutally
            %be different yTrail and yLead from the y1n and y2n
            %y2nTemp(j) = tempTrailR(2);
            zTrailR(end) = tempTrailR(3);
            xTrailL(end) = tempTrailL(1);
            %y2nTemp(j) = tempTrailL(2);
            zTrailL(end) = tempTrailL(3);
            xLeadR(end) = tempLeadR(1);
            %y1nTemp(j) = tempLeadR(2);
            zLeadR(end) = tempLeadR(3);
            xLeadL(end) = tempLeadL(1);
            %y2nTemp(j) = tempLeadL(2);
            zLeadL(end) = tempLeadL(3);
            xm(end) = tempM(1);
            ym(end) = tempM(2);
            zm(end) = tempM(3);
        end
    end
    yTrailR = [yTrailR;wing.Y(1:end-1)+wing.start(2)];
    yTrailL = [yTrailL;wing.Y(2:end)+wing.start(2)];
    yLeadR = [yLeadR;wing.Y(1:end-1)+wing.start(2)];
    yLeadL = [yLeadL;wing.Y(2:end)+wing.start(2)];
    y1n = [y1n;y1nTemp];
    y2n = [y2n;y2nTemp];
end
end