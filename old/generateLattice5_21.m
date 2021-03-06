%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lattice] = generateLattice5_21(SegNum,Z,Phiz,X,Phix,Theta,PhiTheta,cord,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generateLattice5_6: Function for Dynamic TORNADO						
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

%Temp Value for testing
bh = 0.103378;

%Intializing lattice
lattice.COLLOC=[];
lattice.VORTEX=[];
lattice.N=[];
lattice.XYZ=[];

%Initalize Position arrays
x1n = [];
x2n = [];
xTrail = [];
xLead = [];
yTrail = [];
yLead = [];
zTrail = [];
zLead = [];
zTrailn = [];
zLeadn = [];
y1n = [];
y2n = [];
z1n = [];
z2n = [];
x1n = [];
x2n = [];
ym = [];
DiHiAng = [];
xmb = [];
ymb = [];
zmb = [];

xTrailc = [];
xLeadc = [];
zTrailc = [];
zLeadc = [];
zTrailcS = [];
zLeadcS = [];
yTrailcS = [];
yLeadcS = [];

% Create wings
for i = 1:cordNum
    x_3quart = cord/cordNum*(i-.25)-centroid(1);
    x_quart = cord/cordNum*(i-.75)-centroid(1);
    z_3quart = 0;
    z_quart = 0;
    z_trail = 0;
    z_lead = 0;
    for j = 1:SegNum
            DiHiAng(j+SegNum*(i-1)) = atan(.5*(Phiz(j)+Phiz(j+1)));
            alpha(j+SegNum*(i-1)) = 0.5*Theta(j)+.125*L(j)*PhiTheta(j)+0.5*Theta(j+1)-0.125*L(j)*PhiTheta(j+1);
            zm(j+SegNum*(i-1)) = 0.5*Z(j)+.125*Phiz(j)+0.5*Z(j+1)-0.125*Phiz(j+1)+x_3quart*sin(alpha(j))-cos(alpha(j))*z_3quart;
            xm(j+SegNum*(i-1)) = 0.5*X(j)+.125*Phix(j)+0.5*X(j+1)-0.125*Phix(j+1)+x_3quart*cos(alpha(j))+sin(alpha(j))*z_3quart;
    end
    ym = [ym,L(1)/2:L(1):sum(L(1:SegNum/2)),sum(L(1:SegNum/2))+L(SegNum/2)/2+bw:L(1):sum(L)+bw];
    x1n = [x1n;X(1:end-1)+x_quart*cos(Theta(1:end-1))-sin(Theta(1:end-1))*z_quart];
    x2n = [x2n;X(2:end)+x_quart*cos(Theta(2:end))-sin(Theta(2:end))*z_quart];
    xTrail = [xTrail;((cord/cordNum*i-centroid(1))*cos(Theta(2:end))-sin(Theta(2:end))*z_trail)];
    xLead = [xLead;((cord/cordNum*(i-1)-centroid(1))*cos(Theta(2:end))-sin(Theta(2:end))*z_lead)];
    zTrail = [zTrail;zeros(size(Theta(2:end)))];
    zLead = [zLead;zeros(size(Theta(2:end)))];
    y1n = [y1n;0;cumsum(L(1:end/2-1));cumsum(L(end/2:end-1))+sum(L(1:end/2-1))+bw];
    y2n = [y2n;cumsum(L(1:end/2));cumsum(L(end/2+1:end))+sum(L(1:end/2))+bw];
    z1n = [z1n;Z(1:end-1)+sin(Theta(1:end-1))*x_quart+cos(Theta(1:end-1))*z_quart];
    z2n = [z2n;Z(1:end-1)+sin(Theta(2:end))*x_quart+cos(Theta(2:end))*z_quart];
end

%Determine start of nose
if noseH >= bc/cordNum
    noseCord = round(noseH*cordNum/bc);
else
    noseCord = 1;
end

ymb = ym';
xmb = xm';
zmb = [zm';(zb+bh/2)*ones((cordNum-noseCord)*numSpanB,1);(zb-bh/2)*ones((cordNum-noseCord)*numSpanB,1)];
DiHiAng = DiHiAng';

%Create top of the body
for i = noseCord+1:cordNum
    ymb = [ymb;[bw/(2*numSpanB):bw/numSpanB:bw]'+sum(L)/2];
    x_quart = bc/cordNum*(i-.75)-centroid(1)-noseH;
    x_3quart = bc/cordNum*(i-.25)-centroid(1)-noseH;
    xmb = [xmb;x_3quart*ones(numSpanB,1)];
    x1n = [x1n;x_quart*ones(numSpanB,1)];
    x2n = [x2n;x_quart*ones(numSpanB,1)];
    xTrail = [xTrail;(bc/cordNum*i-centroid(1)-noseH)*ones(numSpanB,1)];
    xLead = [xLead;(bc/cordNum*(i-1)-centroid(1)-noseH)*ones(numSpanB,1)];
    y1n = [y1n;[sum(L)/2:bw/numSpanB:sum(L)/2+bw*(1-1/numSpanB)]'];
    y2n = [y2n;[sum(L)/2+bw/numSpanB:bw/numSpanB:sum(L)/2+bw]'];
    z1n = [z1n;(zb+bh/2)*ones(numSpanB,1)];
    z2n = [z2n;(zb+bh/2)*ones(numSpanB,1)];
    zTrail = [zTrail;(zb+bh/2)*ones(numSpanB,1)];
    zLead = [zLead;(zb+bh/2)*ones(numSpanB,1)];
    DiHiAng = [DiHiAng;zeros(numSpanB,1)];
end

%Create Bottom of the body
for i = noseCord+1:cordNum
    ymb = [ymb;[bw/(2*numSpanB):bw/numSpanB:bw]'+sum(L)/2];
    x_quart = bc/cordNum*(i-.75)-centroid(1)-noseH;
    x_3quart = bc/cordNum*(i-.25)-centroid(1)-noseH;
    xmb = [xmb;x_3quart*ones(numSpanB,1)];
    x1n = [x1n;x_quart*ones(numSpanB,1)];
    x2n = [x2n;x_quart*ones(numSpanB,1)];
    xTrail = [xTrail;(bc/cordNum*i-centroid(1)-noseH)*ones(numSpanB,1)];
    xLead = [xLead;(bc/cordNum*(i-1)-centroid(1)-noseH)*ones(numSpanB,1)];
    y1n = [y1n;[sum(L)/2:bw/numSpanB:sum(L)/2+bw*(1-1/numSpanB)]'];
    y2n = [y2n;[sum(L)/2+bw/numSpanB:bw/numSpanB:sum(L)/2+bw]'];
    z1n = [z1n;(zb-bh/2)*ones(numSpanB,1)];
    z2n = [z2n;(zb-bh/2)*ones(numSpanB,1)];
    zTrail = [zTrail;(zb-bh/2)*ones(numSpanB,1)];
    zLead = [zLead;(zb-bh/2)*ones(numSpanB,1)];
    DiHiAng = [DiHiAng;zeros(numSpanB,1)];
end

%Create Left side of the body
for i = noseCord+1:cordNum
    zmb = [zmb;[bh/(2*numSpanB):bh/numSpanB:bh]'+zb-bh/2];
    ymb = [ymb;ones(numSpanB,1)*sum(L)/2];
    x_quart = bc/cordNum*(i-.75)-centroid(1)-noseH;
    x_3quart = bc/cordNum*(i-.25)-centroid(1)-noseH;
    xmb = [xmb;x_3quart*ones(numSpanB,1)];
    x1n = [x1n;x_quart*ones(numSpanB,1)];
    x2n = [x2n;x_quart*ones(numSpanB,1)];
    xTrail = [xTrail;(bc/cordNum*i-centroid(1)-noseH)*ones(numSpanB,1)];
    xLead = [xLead;(bc/cordNum*(i-1)-centroid(1)-noseH)*ones(numSpanB,1)];
    y1n = [y1n;ones(numSpanB,1)*sum(L)/2];
    y2n = [y2n;ones(numSpanB,1)*sum(L)/2];
    z1n = [z1n;[-bh/2:bh/numSpanB:bh/2*(1-1/numSpanB)]'+zb];
    z2n = [z2n;[bh/numSpanB-bh/2:bh/numSpanB:bh/2]'+zb];
    zTrail = [zTrail;[-bh/2:bh/numSpanB:bh/2*(1-1/numSpanB)]'+zb];
    zLead = [zLead;[bh/numSpanB-bh/2:bh/numSpanB:bh/2]'+zb];
    DiHiAng = [DiHiAng;0*ones(numSpanB,1)];
end

%Create Right side of the body
for i = noseCord+1:cordNum
    zmb = [zmb;[bh/(2*numSpanB):bh/numSpanB:bh]'+zb-bh/2];
    ymb = [ymb;ones(numSpanB,1)*sum(L)/2+bw];
    x_quart = bc/cordNum*(i-.75)-centroid(1)-noseH;
    x_3quart = bc/cordNum*(i-.25)-centroid(1)-noseH;
    xmb = [xmb;x_3quart*ones(numSpanB,1)];
    x1n = [x1n;x_quart*ones(numSpanB,1)];
    x2n = [x2n;x_quart*ones(numSpanB,1)];
    xTrail = [xTrail;(bc/cordNum*i-centroid(1)-noseH)*ones(numSpanB,1)];
    xLead = [xLead;(bc/cordNum*(i-1)-centroid(1)-noseH)*ones(numSpanB,1)];
    y1n = [y1n;ones(numSpanB,1)*(sum(L)/2+bw)];
    y2n = [y2n;ones(numSpanB,1)*(sum(L)/2+bw)];
    z1n = [z1n;[-bh/2:bh/numSpanB:bh/2*(1-1/numSpanB)]'+zb];
    z2n = [z2n;[bh/numSpanB-bh/2:bh/numSpanB:bh/2]'+zb];
    zTrail = [zTrail;[-bh/2:bh/numSpanB:bh/2*(1-1/numSpanB)]'+zb];
    zLead = [zLead;[bh/numSpanB-bh/2:bh/numSpanB:bh/2]'+zb];
    DiHiAng = [DiHiAng;pi*ones(numSpanB,1)];
end

temp = [];%change variable name to something meaningful
% Create end cap
% for i = 1:numSpanB
%     zmb = [zmb;[.75*bh/numSpanB:bh/numSpanB:bh/2]'+zb-bh/2];
%     zmb = [zmb;-[.75*bh/numSpanB:bh/numSpanB:bh/2]'+zb+bh/2];
%     ymb = [ymb;bw/numSpanB*(i-.5)*ones(numSpanB,1)+sum(L)/2];
%     xmb = [xmb;(bc-centroid(1)-noseH)*ones(numSpanB,1)];
%     x1n = [x1n;(bc-centroid(1)-noseH)*ones(numSpanB,1)];
%     x2n = [x2n;(bc-centroid(1)-noseH)*ones(numSpanB,1)];
%     xTrail = [xTrail;(bc-centroid(1)-noseH)*ones(numSpanB,1)];
%     xLead = [xLead;(bc-centroid(1)-noseH)*ones(numSpanB,1)];
%     y1n = [y1n;sum(L)/2+bw/numSpanB*(i-1)*ones(numSpanB,1)];
%     y2n = [y2n;sum(L)/2+bw/numSpanB*i*ones(numSpanB,1)];
%     z1n = [z1n;[.25*bh/numSpanB:bh/numSpanB:bh*.5]'+zb-bh/2];
%     z1n = [z1n;-[.25*bh/numSpanB:bh/numSpanB:bh*.5]'+zb+bh/2];
%     z2n = [z2n;[.25*bh/numSpanB:bh/numSpanB:bh*.5]'+zb-bh/2];
%     z2n = [z2n;-[.25*bh/numSpanB:bh/numSpanB:bh*.5]'+zb+bh/2];
%     zTrailn = [zTrailn;[-bh/2:bh/numSpanB:bh/2*(1-1/numSpanB)]'+zb];
%     zLeadn = [zLeadn;[bh/numSpanB-bh/2:bh/numSpanB:bh/2]'+zb];
%     DiHiAng = [DiHiAng;0*ones(numSpanB,1)];
%     temp = [temp,pi*ones(1,numSpanB/2),0*ones(1,numSpanB/2)];
% end

%Create nose flat front
% for i = 1:numSpanB
%     zmb = [zmb;[.75*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb-bh/2];%Simplify equations later
%     zmb = [zmb;-[.75*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb+bh/2];
%     ymb = [ymb;bw/numSpanB*(i-.5)*ones(numSpanB/2,1)+sum(L)/2];
%     xmb = [xmb;(-centroid(1)-noseH)*ones(numSpanB/2,1)];
%     x1n = [x1n;(-centroid(1)-noseH)*ones(numSpanB/2,1)];
%     x2n = [x2n;(-centroid(1)-noseH)*ones(numSpanB/2,1)];
%     xTrail = [xTrail;(-centroid(1)-noseH)*ones(numSpanB/2,1)];
%     xLead = [xLead;(-centroid(1)-noseH)*ones(numSpanB/2,1)];
%     y1n = [y1n;sum(L)/2+bw/numSpanB*(i-1)*ones(numSpanB/2,1)];
%     y2n = [y2n;sum(L)/2+bw/numSpanB*i*ones(numSpanB/2,1)];
%     z1n = [z1n;[.25*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb-bh/2];
%     z1n = [z1n;-[.25*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb+bh/2];
%     z2n = [z2n;[.25*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb-bh/2];
%     z2n = [z2n;-[.25*bh/numSpanB+(numSpanB/4-1)*bh/numSpanB:bh/numSpanB:bh*(.5-1/numSpanB)]'+zb+bh/2];
%     zTrailn = [zTrailn;[-bh/4:bh/numSpanB:bh/4*(1-1/numSpanB)]'+zb];
%     zLeadn = [zLeadn;[bh/numSpanB-bh/4:bh/numSpanB:bh/4]'+zb];
%     DiHiAng = [DiHiAng;0*ones(numSpanB,1)];
%     temp = [temp,0*ones(1,numSpanB/4),pi*ones(1,numSpanB/4)];
% end

%Create nose "cone"
SNose = [];
%Create opperational triangle between flat nose and body in z direction
hyp = sqrt((noseCord*bc/cordNum)^2+(bh/2)^2);
sinn = -bh/(2*hyp);
cosn = noseCord*bc/(cordNum*hyp);
slopeR = [cosn,0,sinn;0,1,0;-sinn,0,cosn];

ymbNose = bw/numSpanB*[.5:numSpanB-.5]+sum(L)/2;
y1nNose = bw/numSpanB*[0:numSpanB-1]+sum(L)/2;
y2nNose = bw/numSpanB*[1:numSpanB]+sum(L)/2;

xmbNose = .75*hyp/noseCord:hyp/noseCord:hyp;
x1nNose = .25*hyp/noseCord:hyp/noseCord:hyp;
x2nNose = .25*hyp/noseCord:hyp/noseCord:hyp;
xLeadNose = 0:hyp/noseCord:hyp-hyp/noseCord;
xTrailNose = hyp/noseCord:hyp/noseCord:hyp;

yLead = y1n-(sum(L)+bw)/2;
yTrail = y2n-(sum(L)+bw)/2;

for i = 0:pi:pi
    rotX = [1,0,0;0,cos(i),-sin(i);0,sin(i),cos(i)];
    
    noseRelmb = slopeR*[combvec(xmbNose,ymbNose);zeros(1,length(ymbNose)*length(xmbNose))];
    noseRel1n = slopeR*[combvec(x1nNose,y1nNose);zeros(1,length(y1nNose)*length(x1nNose))];
    noseRel2n = slopeR*[combvec(x2nNose,y2nNose);zeros(1,length(y2nNose)*length(x2nNose))];
    noseReLeadL = slopeR*[combvec(xLeadNose,y1nNose);zeros(1,length(y1nNose)*length(xLeadNose))];
    noseReLeadR = slopeR*[combvec(xLeadNose,y2nNose);zeros(1,length(y2nNose)*length(xLeadNose))];
    noseReTrailL = slopeR*[combvec(xTrailNose,y1nNose);zeros(1,length(y1nNose)*length(xTrailNose))];
    noseReTrailR = slopeR*[combvec(xTrailNose,y2nNose);zeros(1,length(y2nNose)*length(xTrailNose))];
    
    noseRelmb = rotX*[noseRelmb(1,:);noseRelmb(2,:)-(sum(L)+bw)/2;noseRelmb(3,:)];
    noseRel1n = rotX*[noseRel1n(1,:);noseRel1n(2,:)-(sum(L)+bw)/2;noseRel1n(3,:)];
    noseRel2n = rotX*[noseRel2n(1,:);noseRel2n(2,:)-(sum(L)+bw)/2;noseRel2n(3,:)];
    noseReLeadL = rotX*[noseReLeadL(1,:);noseReLeadL(2,:)-(sum(L)+bw)/2;noseReLeadL(3,:)];
    noseReLeadR = rotX*[noseReLeadR(1,:);noseReLeadR(2,:)-(sum(L)+bw)/2;noseReLeadR(3,:)];
    noseReTrailL = rotX*[noseReTrailL(1,:);noseReTrailL(2,:)-(sum(L)+bw)/2;noseReTrailL(3,:)];
    noseReTrailR = rotX*[noseReTrailR(1,:);noseReTrailR(2,:)-(sum(L)+bw)/2;noseReTrailR(3,:)];
    
    xmb = [xmb;noseRelmb(1,:)'-noseH-centroid(1)];
    ymb = [ymb;noseRelmb(2,:)'+(sum(L)+bw)/2];
    zmb = [zmb;noseRelmb(3,:)'+zb];

    x1n = [x1n;noseRel1n(1,:)'-noseH-centroid(1)];
    y1n = [y1n;noseRel1n(2,:)'+(sum(L)+bw)/2];
    z1n = [z1n;noseRel1n(3,:)'+zb];

    x2n = [x2n;noseRel2n(1,:)'-noseH-centroid(1)];
    y2n = [y2n;noseRel2n(2,:)'+(sum(L)+bw)/2];
    z2n = [z2n;noseRel2n(3,:)'+zb];

    xLead = [xLead;noseReLeadL(1,:)'-noseH-centroid(1)];
    yLead = [yLead;noseRel1n(2,:)'];
    zLeadc = [zLeadc;noseReLeadL(3,:)'+zb];
    
    xTrail = [xTrail;noseReTrailL(1,:)'-noseH-centroid(1)];
    yTrail = [yTrail;noseRel2n(2,:)'];
    zTrailc = [zTrailc;noseReTrailL(3,:)'+zb];

    SNose = [SNose,zeros(1,length(y1nNose)*length(xTrailNose))];
    DiHiAng = [DiHiAng;0*ones(length(y1nNose)*length(xTrailNose),1)];
end

hyp = sqrt((noseCord*bc/cordNum)^2+(bw/2)^2);
sinn = bw/(2*hyp);
cosn = noseCord*bc/(cordNum*hyp);
cost = noseCord*bc/(cordNum*sqrt((noseCord*bc/cordNum)^2+(bh/2)^2));
%slopeR = [cosn,-sinn,0;sinn,cosn,0;0,0,1];
slopeR = eye(3);

zmbNose = bh/numSpanB*[numSpanB/4+.5:3*numSpanB/4-.5]-bh/2;
z1nNose = bh/numSpanB*[numSpanB/4:3*numSpanB/4-1]-bh/2;
z2nNose = bh/numSpanB*[numSpanB/4+1:3*numSpanB/4]-bh/2;

% for i = 0:pi:pi
%     for j = 1:noseCord
%         zmbNose = bh/numSpanB*[numSpanB/(2*noseCord)*j+.5:numSpanB-numSpanB/(2*noseCord)*j-.5]-bh/2;
%         z1nNose = bh/numSpanB*[numSpanB/(2*noseCord)*j:numSpanB-1-numSpanB/(2*noseCord)*j]-bh/2;
%         z2nNose = bh/numSpanB*[numSpanB/(2*noseCord)*j+1:numSpanB-numSpanB/(2*noseCord)*j]-bh/2;
%         
%         xmbNose = noseCord*bc/cordNum-j*bc/cordNum+.75*bc/cordNum;
%         x1nNose = noseCord*bc/cordNum-j*bc/cordNum+.25*bc/cordNum;
%         x2nNose = noseCord*bc/cordNum-j*bc/cordNum+.25*bc/cordNum;
%         xLeadNose = noseCord*bc/cordNum-j*bc/cordNum;
%         xTrailNose = noseCord*bc/cordNum-j*bc/cordNum+bc/cordNum;
% 
%         rotX = [1,0,0;0,cos(i),-sin(i);0,sin(i),cos(i)];
% 
%         tempComb = combvec(xmbNose,zmbNose);
%         noseRelmb = slopeR*[tempComb(1,:);zeros(1,length(zmbNose)*length(xmbNose));tempComb(2,:);];
%         tempComb = combvec(x1nNose,z1nNose);
%         noseRel1n = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%         tempComb = combvec(x2nNose,z2nNose);
%         noseRel2n = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%         tempComb = combvec(xLeadNose,z1nNose);
%         noseReLeadL = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%         tempComb = combvec(xLeadNose,z2nNose);
%         noseReLeadR = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%         tempComb = combvec(xTrailNose,z1nNose);
%         noseReTrailL = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%         tempComb = combvec(xTrailNose,z2nNose);
%         noseReTrailR = slopeR*[tempComb(1,:);zeros(1,length(z1nNose)*length(x1nNose));tempComb(2,:);];
%     
%         noseRelmb = rotX*[noseRelmb(1,:);noseRelmb(2,:)+bw/2;noseRelmb(3,:)];
%         noseRel1n = rotX*[noseRel1n(1,:);noseRel1n(2,:)+bw/2;noseRel1n(3,:)];
%         noseRel2n = rotX*[noseRel2n(1,:);noseRel2n(2,:)+bw/2;noseRel2n(3,:)];
%         noseReLeadL = rotX*[noseReLeadL(1,:);noseReLeadL(2,:)+bw/2;noseReLeadL(3,:)];
%         noseReLeadR = rotX*[noseReLeadR(1,:);noseReLeadR(2,:)+bw/2;noseReLeadR(3,:)];
%         noseReTrailL = rotX*[noseReTrailL(1,:);noseReTrailL(2,:)+bw/2;noseReTrailL(3,:)];
%         noseReTrailR = rotX*[noseReTrailR(1,:);noseReTrailR(2,:)+bw/2;noseReTrailR(3,:)];
% 
%         xmb = [xmb;noseRelmb(1,:)'-noseH-centroid(1)];
%         ymb = [ymb;noseRelmb(2,:)'+(sum(L)+bw)/2];
%         zmb = [zmb;noseRelmb(3,:)'+zb];
% 
%         x1n = [x1n;noseRel1n(1,:)'-noseH-centroid(1)];
%         y1n = [y1n;noseRel1n(2,:)'+(sum(L)+bw)/2];
%         z1n = [z1n;noseRel1n(3,:)'+zb];
% 
%         x2n = [x2n;noseRel2n(1,:)'-noseH-centroid(1)];
%         y2n = [y2n;noseRel2n(2,:)'+(sum(L)+bw)/2];
%         z2n = [z2n;noseRel2n(3,:)'+zb];
% 
%         xLead = [xLead;noseReLeadL(1,:)'-noseH-centroid(1)];
%         yLeadcS = [yLeadcS;noseReLeadL(2,:)'];
%         zLeadcS = [zLeadcS;noseReLeadL(3,:)'+zb];
% 
%         xTrail = [xTrail;noseReTrailL(1,:)'-noseH-centroid(1)];
%         yTrailcS = [yTrailcS;noseReTrailL(2,:)'];
%         zTrailcS = [zTrailcS;noseReTrailR(3,:)'+zb];
% 
%         SNose = [SNose,pi*ones(1,length(z1nNose)*length(x1nNose))];
%         DiHiAng = [DiHiAng;0*ones(length(z1nNose)*length(x1nNose),1)];
%     end
% end

ymb = ymb-(sum(L)+bw)/2;
y1n = y1n-(sum(L)+bw)/2;
y2n = y2n-(sum(L)+bw)/2;
%alpha = alpha+state.alpha;

lattice.COLLOC = [xmb,ymb,zmb];

lattice.VORTEX(:,:,1) = [cord*ones(size(x1n)),x1n,x1n,cord*ones(size(x1n))];
lattice.VORTEX(:,:,2) = [y1n,y1n,y2n,y2n];
lattice.VORTEX(:,:,3) = [z1n,z1n,z2n,z2n];

lattice.XYZ(:,:,1) = [xLead,xLead,xTrail,xTrail,xLead];
lattice.XYZ(:,:,2) = [[yLead;yLeadcS],[yTrail;yLeadcS],[yTrail;yTrailcS],[yLead;yTrailcS],[yLead;yLeadcS]];
lattice.XYZ(:,:,3) = [[zLead;zLeadn;zLeadc;zLeadcS],[zTrail;zLeadn;zLeadc;zLeadcS;],[zTrail;zTrailn;zTrailc;zLeadcS;],[zLead;zTrailn;zTrailc;zTrailcS],[zLead;zLeadn;zLeadc;zTrailcS]];

S = [zeros(1,cordNum*SegNum),0*ones(1,4*(cordNum-noseCord)*numSpanB),temp,SNose];

lattice.N=normals4(lattice.COLLOC,lattice.VORTEX,S,DiHiAng);

[ref]=setRef(L,bw,cord,bc,centroid,noseH);
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
function [normal]=normals4(colloc,vortex,C_Slope,DiHi)
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
   DiHiAngle = DiHi(t);
   
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
         R2=trot3(r0,R,-alpha,DiHiAngle);		%rotating wha trot
         N=[N;R2']; 
end

normal=N;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[p2]=trot3(hinge,p,alpha,DiHi)
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
RZDiHi=[[1 0 0];[0 cos(DiHi) -sin(DiHi)];[0 sin(DiHi) cos(DiHi)]];
RYMT=[[cost 0 -sint];[0 1 0];[sint 0 cost]];
RZMF=[[cosf sinf 0];[-sinf cosf 0];[0 0 1]];

P=RZDiHi*RZF*RYT*RZA*RYMT*RZMF;
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

infx=infdist*cos(state.alpha)*cos(state.betha);
infy=-infdist*sin(state.betha);
infz=infdist*sin(state.alpha)*cos(state.betha);

for t=1:a
   for s=1:2
   	x=infx+lattice.VORTEX(t,c(s),1);
      y=infy+lattice.VORTEX(t,c(s),2);
      z=infz+lattice.VORTEX(t,c(s),3);
      
      psi=state.P/state.AS*x;
  	  theta=state.Q/state.AS*x;
   	  fi=state.R/state.AS*x;
      
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