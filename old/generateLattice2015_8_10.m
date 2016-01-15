%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lattice] = generateLattice2015_8_10(geo,bw,bc,numSpanB,numCordB,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generateLattice6_18_2015: Function for Dynamic TORNADO						
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
y1n = [];
y2n = [];
z1n = [];
z2n = [];
ym = [];
xm = [];
zm = [];

%Iterate through all the wings 
for k = 1:geo.nwing
    %Go through each cord section
    for i = 1:geo.Wings(k).wing.cordNum
        SegCord = [geo.Wings(k).wing.meanCord(i+1,1)-geo.Wings(k).wing.meanCord(i,1);geo.Wings(k).wing.meanCord(i+1,2)-geo.Wings(k).wing.meanCord(i,2)];
        quartCord = .25*SegCord+geo.Wings(k).wing.meanCord(i,:)';
        threeQuartCord = .75*SegCord+geo.Wings(k).wing.meanCord(i,:)';
        trail = geo.Wings(k).wing.meanCord(i+1,:)';
        lead = geo.Wings(k).wing.meanCord(i,:)';
        
        for j = 1:geo.Wings(k).wing.SegNum
            DiHiAng(j+geo.Wings(k).wing.SegNum*(i-1)) = atan(.5*(geo.Wings(k).wing.PhiZ(j+1)+geo.Wings(k).wing.PhiZ(j)));
            alpha(j+geo.Wings(k).wing.SegNum*(i-1)) = 0.5*geo.Wings(k).wing.Theta(j+1)+0.5*geo.Wings(k).wing.Theta(j);
            if j < (geo.Wings(k).wing.SegNum-1)/2+1
                zm = [zm;0.5*geo.Wings(k).wing.Z(j+1)+.125*geo.Wings(k).wing.PhiZ(j+1)+0.5*geo.Wings(k).wing.Z(j)-0.125*geo.Wings(k).wing.PhiZ(j)+[sin(alpha(j)),cos(alpha(j))]*threeQuartCord+geo.Wings(k).wing.start(3)];
                xm = [xm;0.5*geo.Wings(k).wing.X(j+1)+.125*geo.Wings(k).wing.PhiX(j+1)+0.5*geo.Wings(k).wing.X(j)-0.125*geo.Wings(k).wing.PhiX(j)+[cos(alpha(j)),-sin(alpha(j))]*threeQuartCord+geo.Wings(k).wing.start(1)];
            else
                zm = [zm;0.5*geo.Wings(k).wing.Z(j)+.125*geo.Wings(k).wing.PhiZ(j)+0.5*geo.Wings(k).wing.Z(j+1)-0.125*geo.Wings(k).wing.PhiZ(j+1)+[sin(alpha(j)),cos(alpha(j))]*threeQuartCord+geo.Wings(k).wing.start(3)];
                xm = [xm;0.5*geo.Wings(k).wing.X(j)+.125*geo.Wings(k).wing.PhiX(j)+0.5*geo.Wings(k).wing.X(j+1)-0.125*geo.Wings(k).wing.PhiX(j+1)+[cos(alpha(j)),-sin(alpha(j))]*threeQuartCord+geo.Wings(k).wing.start(1)];
            end
        end
        ym = [ym;-flipud([geo.Wings(k).wing.L(1)/2:geo.Wings(k).wing.L(1):sum(geo.Wings(k).wing.L)/2]')-bw/2;[geo.Wings(k).wing.L(1)/2:geo.Wings(k).wing.L(1):sum(geo.Wings(k).wing.L)/2]'+bw/2+geo.Wings(k).wing.start(2)];
        %
        x1n = [x1n;geo.Wings(k).wing.X(1:end-1)+[cos(geo.Wings(k).wing.Theta(1:end-1)),-sin(geo.Wings(k).wing.Theta(1:end-1))]*quartCord+geo.Wings(k).wing.start(1)];
        x2n = [x2n;geo.Wings(k).wing.X(2:end)+[cos(geo.Wings(k).wing.Theta(2:end)),-sin(geo.Wings(k).wing.Theta(2:end))]*quartCord+geo.Wings(k).wing.start(1)];
        xTrailR = [xTrailR;[cos(geo.Wings(k).wing.Theta(1:end-1)),-sin(geo.Wings(k).wing.Theta(1:end-1))]*trail+geo.Wings(k).wing.start(1)];
        xTrailL = [xTrailL;[cos(geo.Wings(k).wing.Theta(2:end)),-sin(geo.Wings(k).wing.Theta(2:end))]*trail+geo.Wings(k).wing.start(1)];
        xLeadR = [xLeadR;[cos(geo.Wings(k).wing.Theta(1:end-1)),-sin(geo.Wings(k).wing.Theta(1:end-1))]*lead+geo.Wings(k).wing.start(1)];
        xLeadL = [xLeadL;[cos(geo.Wings(k).wing.Theta(2:end)),-sin(geo.Wings(k).wing.Theta(2:end))]*lead+geo.Wings(k).wing.start(1)];
        %
        y1n = [y1n;[-flipud([cumsum(geo.Wings(k).wing.L(length(geo.Wings(k).wing.L)/2:end-1))])-bw/2;[0;cumsum(geo.Wings(k).wing.L(length(geo.Wings(k).wing.L)/2:end-2))]+bw/2+geo.Wings(k).wing.start(2)]];
        y2n = [y2n;[-flipud([0;cumsum(geo.Wings(k).wing.L(length(geo.Wings(k).wing.L)/2:end-2))])-bw/2;[cumsum(geo.Wings(k).wing.L(length(geo.Wings(k).wing.L)/2:end-1))]+bw/2+geo.Wings(k).wing.start(2)]];
        %
        z1n = [z1n;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(1:end-1)),cos(geo.Wings(k).wing.Theta(1:end-1))]*quartCord+geo.Wings(k).wing.start(3)];
        z2n = [z2n;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(2:end)),cos(geo.Wings(k).wing.Theta(2:end))]*quartCord+geo.Wings(k).wing.start(3)];
        zTrailR = [zTrailR;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(1:end-1)),cos(geo.Wings(k).wing.Theta(1:end-1))]*trail+geo.Wings(k).wing.start(3)];
        zTrailL = [zTrailL;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(2:end)),cos(geo.Wings(k).wing.Theta(2:end))]*trail+geo.Wings(k).wing.start(3)];
        zLeadR = [zLeadR;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(1:end-1)),cos(geo.Wings(k).wing.Theta(1:end-1))]*lead+geo.Wings(k).wing.start(3)];
        zLeadL = [zLeadL;geo.Wings(k).wing.Z(1:end-1)+[sin(geo.Wings(k).wing.Theta(2:end)),cos(geo.Wings(k).wing.Theta(2:end))]*lead+geo.Wings(k).wing.start(3)];
    end
end
ymb = ym;
xmb = xm;
zmb = [zm;zeros(numCordB*numSpanB,1)];
DiHiAng = DiHiAng';

if numSpanB ~= 0 || bw ~=0
    for i = 1:numCordB
        ymb = [ymb;[-flipud([bw/(2*numSpanB):bw/numSpanB:bw/2]');[bw/(2*numSpanB):bw/numSpanB:bw/2]']];
        x_quart = bc/numCordB*(i-.75);
        x_3quart = bc/numCordB*(i-.25);
        xmb = [xmb;x_3quart*ones(numSpanB,1)];
        x1n = [x1n;x_quart*ones(numSpanB,1)];
        x2n = [x2n;x_quart*ones(numSpanB,1)];
        xTrailR = [xTrailR;(bc/numCordB*i)*ones(numSpanB,1)];
        xTrailL = [xTrailL;(bc/numCordB*i)*ones(numSpanB,1)];
        xLeadR = [xLeadR;(bc/numCordB*(i-1))*ones(numSpanB,1)];
        xLeadL = [xLeadL;(bc/numCordB*(i-1))*ones(numSpanB,1)];
        y1n = [y1n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2]');[0:bw/numSpanB:bw/2-bw/numSpanB]']];
        y2n = [y2n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2-bw/numSpanB]');[0:bw/numSpanB:bw/2]']];
        z1n = [z1n;zeros(numSpanB,1)];
        z2n = [z2n;zeros(numSpanB,1)];
        zTrailR = [zTrailR;zeros(numSpanB,1)];
        zLeadR = [zLeadR;zeros(numSpanB,1)];
        zTrailL = [zTrailL;zeros(numSpanB,1)];
        zLeadL = [zLeadL;zeros(numSpanB,1)];
        DiHiAng = [DiHiAng;zeros(numSpanB,1)];
    end
end

lattice.COLLOC = [xmb,ymb,zmb];

lattice.VORTEX(:,:,1) = [xTrailR,x1n,x1n,xTrailR];
lattice.VORTEX(:,:,2) = [y1n,y1n,y2n,y2n];
lattice.VORTEX(:,:,3) = [z1n,z1n,z2n,z2n];

lattice.XYZ(:,:,1) = [xLeadR,xLeadL,xTrailL,xTrailR,xLeadR];
lattice.XYZ(:,:,2) = [y1n,y2n,y2n,y1n,y1n];
lattice.XYZ(:,:,3) = [zLeadR,zLeadL,zTrailL,zTrailR,zLeadR];

S = zeros(1,geo.nx*geo.ny+numCordB*numSpanB);
lattice.N=normals4(lattice.COLLOC,lattice.VORTEX,S);

[ref]=setRef6_15(geo.Wings(1).wing.span,geo.Wings(1).wing.cord,geo.CG);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function [lattice]=wakesetup2(lattice,state,ref)
infdist=config('infinity');
if isempty(infdist)
	infdist=6*ref.b_ref;
end

[a b c]=size(lattice.VORTEX);
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