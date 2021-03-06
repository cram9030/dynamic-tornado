%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lattice] = generateLattice6_24_2015(SegNum,Z,Phiz,X,Phix,Theta,cord,geo,cordNum,L,bw,bc,noseH,zb,numSpanB,state)
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
% Inputs:   SegNum - Number of segments for both wings
%           Z - Z displacement of wing from centroid
%           Phiz - Z displacement slope
%           X - X displacement of wing from centroid
%           Phix - X displacement slope
%           Theta - Twist angle
%           PhiTheta - Twist angle slope
%           cord - Wing Cord Length
%           centroid - Center of mass of airfoil
%           cordNum - Number of lattice partitions in the cord
%           L - Length vector for wings
%           bw - Body Width
%           bc - Body Length
%           noseH - length that the nose cone extends beyond the wing
%           zb - Z offset of the body from the wing
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
%bh = 0.1034;
bh = 0;
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

for i = 1:cordNum
    x_3quart = cord/cordNum*(i-.25);
    x_quart = cord/cordNum*(i-.75);
    z_3quart = -geo.ref_point(2);
    z_quart = -geo.ref_point(2);
    z_trail = -geo.ref_point(2);
    z_lead = -geo.ref_point(2);
    for j = 1:SegNum
        DiHiAng(j+SegNum*(i-1)) = atan(.5*(Phiz(j+1)+Phiz(j)));
        alpha(j+SegNum*(i-1)) = 0.5*Theta(j+1)+0.5*Theta(j);
        if j < (SegNum-1)/2+1
            zm(j+SegNum*(i-1)) = 0.5*Z(j+1)+.125*Phiz(j+1)+0.5*Z(j)-0.125*Phiz(j)+x_3quart*sin(alpha(j))-cos(alpha(j))*z_3quart;
            xm(j+SegNum*(i-1)) = 0.5*X(j+1)+.125*Phix(j+1)+0.5*X(j)-0.125*Phix(j)+x_3quart*cos(alpha(j))+sin(alpha(j))*z_3quart;
        else
            zm(j+SegNum*(i-1)) = 0.5*Z(j)+.125*Phiz(j)+0.5*Z(j+1)-0.125*Phiz(j+1)+x_3quart*sin(alpha(j))-cos(alpha(j))*z_3quart;
            xm(j+SegNum*(i-1)) = 0.5*X(j)+.125*Phix(j)+0.5*X(j+1)-0.125*Phix(j+1)+x_3quart*cos(alpha(j))+sin(alpha(j))*z_3quart;
        end
    end
    ym = [ym;-flipud([L(1)/2:L(1):sum(L)/2]')-bw/2;[L(1)/2:L(1):sum(L)/2]'+bw/2];
    x1n = [x1n;X(1:end-1)+x_quart*cos(Theta(1:end-1))-sin(Theta(1:end-1))*z_quart];
    x2n = [x2n;X(2:end)+x_quart*cos(Theta(2:end))-sin(Theta(2:end))*z_quart];
    xTrailR = [xTrailR;((cord/cordNum*i)*cos(Theta(1:end-1))-sin(Theta(1:end-1))*z_trail)];
    xTrailL = [xTrailL;((cord/cordNum*i)*cos(Theta(2:end))-sin(Theta(2:end))*z_trail)];
    xLeadR = [xLeadR;((cord/cordNum*(i-1))*cos(Theta(1:end-1))-sin(Theta(1:end-1))*z_lead)];
    xLeadL = [xLeadL;((cord/cordNum*(i-1))*cos(Theta(2:end))-sin(Theta(2:end))*z_lead)];
    y1n = [y1n;[-flipud([cumsum(L(length(L)/2:end-1))])-bw/2;[0;cumsum(L(length(L)/2:end-2))]+bw/2]];
    y2n = [y2n;[-flipud([0;cumsum(L(length(L)/2:end-2))])-bw/2;[cumsum(L(length(L)/2:end-1))]+bw/2]];
    z1n = [z1n;Z(1:end-1)+sin(Theta(1:end-1))*x_quart+cos(Theta(1:end-1))*z_quart];
    z2n = [z2n;Z(1:end-1)+sin(Theta(2:end))*x_quart+cos(Theta(2:end))*z_quart];
    zTrailR = [zTrailR;Z(1:end-1)+sin(Theta(1:end-1))*(cord/cordNum*i-geo.ref_point(1))+cos(Theta(1:end-1))*z_trail];
    zTrailL = [zTrailL;Z(1:end-1)+sin(Theta(2:end))*(cord/cordNum*i-geo.ref_point(1))+cos(Theta(2:end))*z_trail];
    zLeadR = [zLeadR;Z(1:end-1)+sin(Theta(1:end-1))*(cord/cordNum*(i-1)-geo.ref_point(1))+cos(Theta(1:end-1))*z_lead];
    zLeadL = [zLeadL;Z(1:end-1)+sin(Theta(2:end))*(cord/cordNum*(i-1)-geo.ref_point(1))+cos(Theta(2:end))*z_lead];
end
ymb = ym;
xmb = xm';
zmb = [zm';(zb-bh/4)*ones(cordNum*numSpanB,1)];
DiHiAng = DiHiAng';

for i = 1:cordNum
    ymb = [ymb;[-flipud([bw/(2*numSpanB):bw/numSpanB:bw/2]');[bw/(2*numSpanB):bw/numSpanB:bw/2]']];
    x_quart = bc/cordNum*(i-.75);
    x_3quart = bc/cordNum*(i-.25);
    xmb = [xmb;x_3quart*ones(numSpanB,1)];
    x1n = [x1n;x_quart*ones(numSpanB,1)];
    x2n = [x2n;x_quart*ones(numSpanB,1)];
    xTrailR = [xTrailR;(bc/cordNum*i)*ones(numSpanB,1)];
    xTrailL = [xTrailL;(bc/cordNum*i)*ones(numSpanB,1)];
    xLeadR = [xLeadR;(bc/cordNum*(i-1))*ones(numSpanB,1)];
    xLeadL = [xLeadL;(bc/cordNum*(i-1))*ones(numSpanB,1)];
    y1n = [y1n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2]');[0:bw/numSpanB:bw/2-bw/numSpanB]']];
    y2n = [y2n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2-bw/numSpanB]');[0:bw/numSpanB:bw/2]']];
    z1n = [z1n;(zb-bh/4)*ones(numSpanB,1)];
    z2n = [z2n;(zb-bh/4)*ones(numSpanB,1)];
    zTrailR = [zTrailR;(zb-bh/4)*ones(numSpanB,1)];
    zLeadR = [zLeadR;(zb-bh/4)*ones(numSpanB,1)];
    zTrailL = [zTrailL;(zb-bh/4)*ones(numSpanB,1)];
    zLeadL = [zLeadL;(zb-bh/4)*ones(numSpanB,1)];
    DiHiAng = [DiHiAng;zeros(numSpanB,1)];
end

% zmb = [zmb;(zb+bh/4)*ones(cordNum*numSpanB,1)];
% 
% for i = 1:cordNum
%     ymb = [ymb;[-flipud([bw/(2*numSpanB):bw/numSpanB:bw/2]');[bw/(2*numSpanB):bw/numSpanB:bw/2]']];
%     x_quart = bc/cordNum*(i-.75)-noseH;
%     x_3quart = bc/cordNum*(i-.25)-noseH;
%     xmb = [xmb;x_3quart*ones(numSpanB,1)];
%     x1n = [x1n;x_quart*ones(numSpanB,1)];
%     x2n = [x2n;x_quart*ones(numSpanB,1)];
%     xTrail = [xTrail;(bc/cordNum*i-noseH)*ones(numSpanB,1)];
%     xLead = [xLead;(bc/cordNum*(i-1)-noseH)*ones(numSpanB,1)];
%     y1n = [y1n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2]');[0:bw/numSpanB:bw/2-bw/numSpanB]']];
%     y2n = [y2n;[-flipud([bw/numSpanB:bw/numSpanB:bw/2-bw/numSpanB]');[0:bw/numSpanB:bw/2]']];
%     z1n = [z1n;(zb+bh/4)*ones(numSpanB,1)];
%     z2n = [z2n;(zb+bh/4)*ones(numSpanB,1)];
%     zTrailR = [zTrailR;(zb+bh/4)*ones(numSpanB,1)];
%     zLeadR = [zLeadR;(zb+bh/4)*ones(numSpanB,1)];
%     zTrailL = [zTrailL;(zb+bh/4)*ones(numSpanB,1)];
%     zLeadL = [zLeadL;(zb+bh/4)*ones(numSpanB,1)];
%     DiHiAng = [DiHiAng;zeros(numSpanB,1)];
% end

lattice.COLLOC = [xmb,ymb,zmb];

lattice.VORTEX(:,:,1) = [cord*ones(size(xTrailR)),x1n,x1n,cord*ones(size(xTrailR))];
lattice.VORTEX(:,:,2) = [y1n,y1n,y2n,y2n];
lattice.VORTEX(:,:,3) = [z1n,z1n,z2n,z2n];

lattice.XYZ(:,:,1) = [xLeadR,xLeadL,xTrailL,xTrailR,xLeadR];
lattice.XYZ(:,:,2) = [y1n,y2n,y2n,y1n,y1n];
lattice.XYZ(:,:,3) = [zLeadR,zLeadL,zTrailL,zTrailR,zLeadR];

S = zeros(1,cordNum*(SegNum+numSpanB));
lattice.N=normals4(lattice.COLLOC,lattice.VORTEX,S);

[ref]=setRef6_15(L,cord,geo.CG);
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
% ref: 	R�de, Westergren, BETA 4th ed,   
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