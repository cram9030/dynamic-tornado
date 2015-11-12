function [results]=coeff_create(results,lattice,state,ref,geo) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient creator: Essential function for TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes aerodynamic coefficients			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautical 
%                               and Vehicle Engineering	
%			copyright 2003											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidary function for TORNADO					
% Called by:	solverloop											
% Calls:			MATLAB standard functions																			
% Loads:																
% Saves: 												
% Input: 			
% Output: forces moments coefficients							
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=0.01;
q=0.5*state.rho*state.U_inf^2;			    %calculating dynamic pressure
										%for coefficient calculation		
[a,b,~]=size(results.F);

npan =  zeros(geo.nwing,1);
for i = 1:geo.nwing
    npan(i) = geo.Wings(i).wing.SegNum*geo.Wings(i).wing.cordNum;
end

normal_force = zeros(a,1);
for s=1:a
	normal_force(s)=results.F(s,:)*lattice.N(s,:)';                                
end                                
panel_area=tarea(lattice.XYZ)';
stat_press=normal_force./panel_area;	%Delta pressure, top/bottom
results.cp=((stat_press)./(q))';

sonic=find(results.cp<fSonicCP6_24_2015(state));
results.sonicpanels=zeros(size(results.cp));
if isempty(sonic)
    %Flow is subsonic
    results.sonicWarning=0;
else
    %Flow is partially supersonic
    %tdisp('Supersonic flow detected')
    results.sonicCP=fSonicCP6_24_2015(state);
    results.sonicWarning=1; 
    results.sonicFraction=size(sonic)/size(lattice.N); %ratio of supersonic panels
    results.sonicpanels(sonic)=1;
end
    


results.CX=results.FORCE(1)/(q*ref.S_ref);
results.CY=results.FORCE(2)/(q*ref.S_ref);
results.CZ=results.FORCE(3)/(q*ref.S_ref);

B2WTransform=[cos(state.betha)*cos(state.alpha_root),        -sin(state.betha),          cos(state.betha)*sin(state.alpha_root) ;...
              cos(state.alpha_root)*sin(state.betha),         cos(state.betha),          sin(state.betha)*sin(state.alpha_root) ;...
                              -sin(state.alpha_root),                        0,                           cos(state.alpha_root)];

lemma=B2WTransform*results.FORCE';

results.D=lemma(1);
results.C=lemma(2);
results.L=lemma(3);

results.CL=results.L/(q*ref.S_ref);
results.CD=results.D/(q*ref.S_ref);
results.CC=results.C/(q*ref.S_ref);

results.Cl=results.MOMENTS(1)/(q*ref.S_ref*ref.b_ref);
results.Cm=results.MOMENTS(2)/(q*ref.S_ref*ref.C_mac);
results.Cn=results.MOMENTS(3)/(q*ref.S_ref*ref.b_ref);



%% ------------ CL per wing computation
%npan=cumsum(sum((geo.nx+geo.fnx).*geo.ny,2).*(geo.symetric+1)'); %Number
%of panels per wing, only one wing right now

index1=1;

for i=1:geo.nwing
    index2=npan(i);
    
    lemma2=B2WTransform*(sum(results.F(index1:index2,:)))';
    
    results.CLwing(i)=lemma2(3)/(q*ref.S_ref);
    results.CDwing(i)=lemma2(1)/(q*ref.S_ref);
    results.CYwing(i)=lemma2(2)/(q*ref.S_ref);
    
    index1=npan(i)+1;
end
%% ----------

%%Differentiation
lemma=B2WTransform*results.dFORCE_alpha';
results.CL_a=lemma(1)/(q*ref.S_ref);
results.CD_a=lemma(2)/(q*ref.S_ref);
results.CC_a=lemma(3)/(q*ref.S_ref);

results.Cl_a=results.dMOMENTS_alpha(1)/(q*ref.S_ref*ref.b_ref);
results.Cm_a=results.dMOMENTS_alpha(2)/(q*ref.S_ref*ref.C_mac);
results.Cn_a=results.dMOMENTS_alpha(3)/(q*ref.S_ref*ref.b_ref);

lemma=B2WTransform*results.dFORCE_beta';
results.CL_b=lemma(1)/(q*ref.S_ref);
results.CD_b=lemma(2)/(q*ref.S_ref);
results.CC_b=lemma(3)/(q*ref.S_ref);

results.Cl_b=results.dMOMENTS_beta(1)/(q*ref.S_ref*ref.b_ref);
results.Cm_b=results.dMOMENTS_beta(2)/(q*ref.S_ref*ref.C_mac);
results.Cn_b=results.dMOMENTS_beta(3)/(q*ref.S_ref*ref.b_ref);

%return   
[results]=spanload6(results,geo,lattice,state);

%[lemma]=fStripforce(geo,results,lattice,state,ref,vCfraction)

end%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results]=spanload6(results,geo,lattice,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	CONFIG: Basic computation function   	%		 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Computes the spanload (force/meter) for 
%  all wings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Tomas Melin, KTH, Department of% 
%	Aeronautics, copyright 2002				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Context: Auxillary function for TORNADO%
%	Called by: TORNADO SOlverloop          %
%	Calls:	None									%
%	Loads:	None									%
%	Generates:	force per meter array 
%     			(ystations X wings)			
%					Ystation array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Revision history post alfa 1.0			%
%  2007-02-14  rho moved to state 
%  2002-05-02
%   input var T (taper) added to get local
%	 chords.
% input var AS (airspeed) added 
%   local chord computation function call added
%

lemma=size(results.F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alertnative solution
for i=1:lemma(1)
    forceMagn(i)=-results.F(i,:)*lattice.N(i,:)'; %Force magnitude (3Dvector -> scalar)
   										          %Aligned with panel normals
                                                  
    B2WTransform=[cos(state.betha)*cos(state.alpha(i)),        -sin(state.betha),          cos(state.betha)*sin(state.alpha(i)) ;...
              cos(state.alpha(i))*sin(state.betha),         cos(state.betha),          sin(state.betha)*sin(state.alpha(i)) ;...
                              -sin(state.alpha(i)),                        0,                           cos(state.alpha(i))];
    lemma4(i,:)=B2WTransform*results.F(i,:)';                                         
	forceLift(i)=lemma4(i,3);                     %Lift on each panel, this is outdata for the
                                                  %viscous correction.
end

A1=((lattice.XYZ(:,1,:)-lattice.XYZ(:,2,:)));
p_span=sqrt(A1(:,:,2).^2+A1(:,:,3).^2); %span of each panel
FPM=forceMagn'./p_span;					%Force per meter on each panel.
LPM=forceLift'./p_span;					%Lift per meter on each panel.

for j = 1:geo.nwing
    cord = zeros(geo.Wings(j).wing.SegNum,1);
    tempForce = zeros(geo.Wings(j).wing.SegNum,1);
    shearforce = zeros(geo.Wings(j).wing.SegNum,1);
    bendingmoment = zeros(geo.Wings(j).wing.SegNum,1);
    ForcePerMeter = zeros(geo.Wings(j).wing.SegNum,1);
    CL_local = zeros(geo.Wings(j).wing.SegNum,1);

    for i = 1:geo.Wings(j).wing.SegNum
        ForcePerMeter(i,1) = sum(FPM(i:geo.Wings(j).wing.SegNum:geo.Wings(j).wing.cordNum*geo.Wings(j).wing.SegNum));
        cord(i) = max(max(lattice.XYZ(i:geo.Wings(j).wing.SegNum:geo.Wings(j).wing.cordNum*geo.Wings(j).wing.SegNum,:,1)))-min(min(lattice.XYZ(i:geo.Wings(j).wing.SegNum:geo.Wings(j).wing.cordNum*geo.Wings(j).wing.SegNum,:,1)));
        tempForce(i) = sum(forceMagn(i:geo.Wings(j).wing.SegNum:geo.Wings(j).wing.cordNum*geo.Wings(j).wing.SegNum));
        CL_local(i,1) = sum(LPM(i:geo.Wings(j).wing.SegNum:geo.Wings(j).wing.cordNum*geo.Wings(j).wing.SegNum));
    end

    temp = [lattice.COLLOC(1:geo.Wings(j).wing.SegNum,2),ForcePerMeter,tempForce,CL_local];
    temp = sortrows(temp);

    shearforce(1:(geo.Wings(j).wing.SegNum)/2) = cumsum(temp(1:(geo.Wings(j).wing.SegNum)/2,3));
    shearforce((geo.Wings(j).wing.SegNum)/2+1:end) = -flipud(cumsum(flipud(temp((geo.Wings(j).wing.SegNum)/2+1:end,3))));


    bendingmoment(1:(geo.Wings(j).wing.SegNum)/2) = cumsum(temp(1:(geo.Wings(j).wing.SegNum)/2,3).*temp(1:(geo.Wings(j).wing.SegNum)/2,1));
    bendingmoment((geo.Wings(j).wing.SegNum)/2+1:end) = flipud(cumsum(-flipud(temp((geo.Wings(j).wing.SegNum)/2+1:end,3)).*flipud(temp((geo.Wings(j).wing.SegNum)/2+1:end,1))));
    
    results.spanResults(j).ystation = temp(:,1);
    results.spanResults(j).ForcePerMeter = temp(:,2);
    results.spanResults(j).ShearForce = shearforce;
    results.spanResults(j).BendingMoment = bendingmoment;
    results.spanResults(j).CL_local = temp(:,3);
end

end%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [panel_area]=tarea(XYZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tarea: Subsidary function for TORNADO					   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the area of each panel								
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautics	%
%				Copyright 2000											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidaty function for TORNADO					
% Called by:	coeff_create
% 
% Calls:			MATLAB 5.2 std fcns								
% Loads:	none
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b c]=size(XYZ);
for i=1:a
   p1=[XYZ(i,1,1) XYZ(i,1,2) XYZ(i,1,3)];	%sets up the vectors 
   p2=[XYZ(i,2,1) XYZ(i,2,2) XYZ(i,2,3)];	%to the corners of the		
   p3=[XYZ(i,3,1) XYZ(i,3,2) XYZ(i,3,3)];	%panel.
   p4=[XYZ(i,4,1) XYZ(i,4,2) XYZ(i,4,3)];
   
   a=p2-p1;	%sets up the edge vectors
   b=p4-p1;
   c=p2-p3;
   d=p4-p3;
   
   ar1=norm(cross(b,a))/2;	%claculates the ctoss product of
   ar2=norm(cross(c,d))/2;	%two diagonal corners
   
 	panel_area(i)=ar1+ar2;	%Sums up the product to make the
end						    %Area
end% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[lc]=fLocal_chord2(geo,lattice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Geometry function 						 	%		 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Computes the Local chord at each collocation 
%  point row.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Tomas Melin, KTH, Department of% 
%	Aeronautics, copyright 2002				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Context: Auxillary function for TORNADO%
%	Called by: TORNADO spanload            %
%	Calls:	None									%
%	Loads:	None									%
%	Generates:	Local chord vector lc, same 
%  order as colloc, N, and the others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[indx1 indx2]=size(geo.b);

lc=[];	%Local chord vector.


panelchords1=sqrt(sum((lattice.XYZ(:,1,:)-lattice.XYZ(:,4,:)).^2,3)); %inboard
panelchords2=sqrt(sum((lattice.XYZ(:,2,:)-lattice.XYZ(:,3,:)).^2,3)); %outboard
panelchords3=(panelchords1+panelchords2)/2; %Chord of each panel, CAUTION 
                                            %this is really camber line
                                            %length, so not really chord
                                            %for very cambered profiles

for i=1:indx1;			%Wing	
   for j=1:indx2;		%Partition
      lemma=[]; %local chord lemma vector.
      chordwisepanels=geo.nx(i,j); %number of panels chordwise on 
                                                %this partition 
      for k=1:geo.ny(i,j)                       %loop over panel strips.
          if geo.ny(i,j)~=0
              lemma=[lemma sum(panelchords3(1:chordwisepanels))];
              panelchords3=panelchords3((chordwisepanels+1):end);
              %size(panelchords3);
          end
      end  
      if geo.symetric(i)==1	%symmetric wings got two sides
         lc=[lc lemma lemma];
         panelchords3=panelchords3((chordwisepanels*geo.ny(i,j)+1):end);
      else
         lc=[lc lemma];
      end
          
   end
end
end%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   