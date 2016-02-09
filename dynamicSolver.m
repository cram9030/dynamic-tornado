function [results]=dynamicSolver(state,geo,lattice,latticei,h,derivative)
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
[rhs,rhsi] = setDynamicBoundary(lattice,state,geo,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving for rhs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

gamma=w2\rhs';
gammai=w2\rhsi;

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



time = cputime;
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

le=(p2-p1);					%Vortex span vector

Lle=sqrt(sum(le.^2,2));
lehat(:,1)=le(:,1)./Lle;
lehat(:,2)=le(:,2)./Lle;
lehat(:,3)=le(:,3)./Lle;


IW(:,1)=DWX*gamma;
IW(:,2)=DWY*gamma;
IW(:,3)=DWZ*gamma;

G(:,1)=gamma.*lehat(:,1);	%Aligning vorticity along panel vortex
G(:,2)=gamma.*lehat(:,2);
G(:,3)=gamma.*lehat(:,3);

V=state.AS;
Wind=([V.*cos(state.alpha).*cos(state.betha) -V.*cos(state.alpha).*sin(state.betha) V.*sin(state.alpha)])-reshape(IW,[length(IW),3]);
for i=1:a
    Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]); %Calculating rotations
end                                   

Wind=Wind+Rot;								%Adding rotations

Fprim=state.rho*cross(Wind,G);			    %Force per unit length
F=Fprim.*(Lle*ones(1,3)); 

results.F=F;
results.FORCE=sum(F,1);						%Total force
M=cross(c3,F);			                 %Moments per panel
results.M=M;
results.MOMENTS=sum(M,1);					%Summing up moments	
results.gamma=gamma;


%Finding derivatives of wind speed, angle of attack and beta via complex
%step function
IWi(:,:,1)=DWX*gammai;
IWi(:,:,2)=DWY*gammai;
IWi(:,:,3)=DWZ*gammai;

G(:,:,1)=gammai.*(lehat(:,1)*ones(1,3));	%Aligning vorticity along panel vortex
G(:,:,2)=gammai.*(lehat(:,2)*ones(1,3));
G(:,:,3)=gammai.*(lehat(:,3)*ones(1,3));

V=state.AS;
Vi=state.AS+1j*h;

if derivative
    WindVi =([Vi.*cos(state.alpha).*cos(state.betha) -Vi.*cos(state.alpha).*sin(state.betha) Vi.*sin(state.alpha)])-reshape(IWi(:,1,:),[length(IWi(:,1,:)),3]);
    WindAlphai = ([V.*cos(state.alpha+1j*h).*cos(state.betha) -V.*cos(state.alpha+1j*h).*sin(state.betha) V.*sin(state.alpha+1j*h)])-reshape(IWi(:,2,:),[length(IWi(:,2,:)),3]);
    WindBetai = ([V.*cos(state.alpha).*cos(state.betha+1j*h) -V.*cos(state.alpha).*sin(state.betha+1j*h) V.*sin(state.alpha)])-reshape(IWi(:,3,:),[length(IWi(:,3,:)),3]);

    results.dF_V = imag(state.rho*cross(WindVi,reshape(G(:,1,:),[length(G(:,1,:)),3])).*(Lle*ones(1,3)))/h;
    results.dFORCE_V = imag(sum(state.rho*cross(WindVi,reshape(G(:,1,:),[length(G(:,1,:)),3])).*(Lle*ones(1,3))))/h;
    results.dM_V = imag(cross(c3,state.rho*cross(WindVi,reshape(G(:,1,:),[length(G(:,1,:)),3])).*(Lle*ones(1,3))))/h;
    results.dMOMENTS_V = imag(sum(cross(c3,state.rho*cross(WindVi,reshape(G(:,1,:),[length(G(:,1,:)),3])).*(Lle*ones(1,3)))))/h;

    results.dF_alpha = imag(state.rho*cross(WindAlphai,reshape(G(:,2,:),[length(G(:,2,:)),3])).*(Lle*ones(1,3)))/h;
    results.dFORCE_alpha = imag(sum(state.rho*cross(WindAlphai,reshape(G(:,2,:),[length(G(:,2,:)),3])).*(Lle*ones(1,3))))/h;
    results.dM_alpha = imag(cross(c3,state.rho*cross(WindAlphai,reshape(G(:,2,:),[length(G(:,2,:)),3])).*(Lle*ones(1,3))))/h;
    results.dMOMENTS_alpha = imag(sum(cross(c3,state.rho*cross(WindAlphai,reshape(G(:,2,:),[length(G(:,2,:)),3])).*(Lle*ones(1,3)))))/h;

    results.NP = pinv([0,results.dFORCE_alpha(3),-results.dFORCE_alpha(2);...
        -results.dFORCE_alpha(3),0,results.dFORCE_alpha(1);...
        results.dFORCE_alpha(2),results.dFORCE_alpha(1),0])*sum(cross(lattice.COLLOC,results.dF_alpha))';

    results.dF_beta = imag(state.rho*cross(WindBetai,reshape(G(:,3,:),[length(G(:,3,:)),3])).*(Lle*ones(1,3)))/h;
    results.dFORCE_beta = imag(sum(state.rho*cross(WindBetai,reshape(G(:,3,:),[length(G(:,3,:)),3])).*(Lle*ones(1,3))))/h;
    results.dM_beta = imag(cross(c3,state.rho*cross(WindBetai,reshape(G(:,3,:),[length(G(:,3,:)),3])).*(Lle*ones(1,3))))/h;
    results.dMOMENTS_beta = imag(sum(cross(c3,state.rho*cross(WindBetai,reshape(G(:,3,:),[length(G(:,3,:)),3])).*(Lle*ones(1,3)))))/h;

    %Calculating flexible derivatives via the complex step method
    for i = 1:length(latticei)
        [w2 void]=fastdw(latticei(i));
        gamma=w2\rhs';

        p1 = [];
        p2 = [];
        p1(:,:)=latticei(i).VORTEX(:,b1,:);		%Calculating panel vortex midpoint	
        p2(:,:)=latticei(i).VORTEX(:,b1+1,:);	%to use as a force locus
        latticei(i).COLLOC(:,:)=(p1+p2)./2;	    % LOCAL control point, vortex midpoint.
        c3=latticei(i).COLLOC-ones(size(latticei(i).COLLOC,1),1)*geo.CG;

        [w3 DW]=fastdw(latticei(i));	                    %Calculating downwash on vorticies

        w4=sum(DW,2);					                %superpositioning aerodynamic influence

        DWX=DW(:,:,1);
        DWY=DW(:,:,2);
        DWZ=DW(:,:,3);

        IW = [];
        IW(:,1)=DWX*gamma;
        IW(:,2)=DWY*gamma;
        IW(:,3)=DWZ*gamma;

        Gi = [];
        Gi(:,1)=gamma.*lehat(:,1);	%Aligning vorticity along panel vortex
        Gi(:,2)=gamma.*lehat(:,2);
        Gi(:,3)=gamma.*lehat(:,3);

        V=state.AS;
        Wind=([V.*cos(state.alpha).*cos(state.betha) -V.*cos(state.alpha).*sin(state.betha) V.*sin(state.alpha)])-reshape(IW,[length(IW),3]);

        Wind=Wind;								%Adding rotations

        Fprim=state.rho*cross(Wind,Gi);			    %Force per unit length
        F=Fprim.*(Lle*ones(1,3));

        dInput(i).dF_control=imag(F)/h;
        dInput(i).dFORCE_control=imag(sum(F,1))/h;						%Total force
        M=cross(c3,F);			                 %Moments per panel
        dInput(i).dM_control=imag(M)/h;
        dInput(i).dMOMENTS_control=imag(sum(M,1))/h;					%Summing up moments	
        dInput(i).dgamma_control=imag(gamma)/h;
    end
    results.dInput = dInput;
end

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
    
    lr1 = [];
    lr2 = [];
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