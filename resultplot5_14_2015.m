function []=resultplot5_14_2015(JID,results,geo,lattice,state,ref);
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
%%% Result plot function
%
% This function plots the results of a tornado run using standard Matlab
% plots. The function reads the output files and displays the JID.
%
% usage:   []=resultplot(VAR)
%
% VAR is the selector for which type of results should be plotted. The
% selector is choosen in postproc, and the list of choises are listed in
% questions(8).
%
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado interface function.
%
% Revision History:
%   Bristol,  2007 06 27:  Addition of new header. TM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings=config('startup');

fname=strcat(JID,'-Cx');
cd(settings.odir)
save(fname,'results','geo','lattice','state','ref')
cd(settings.hdir)

cd(settings.odir)
load(fname)
cd(settings.hdir)   

[x y z]=midpoint(lattice.XYZ);
d=size(z,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)% 	Delta cp plot
rotate3d on 
colormap(hot);
fill3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)',results.cp')
title('Delta cp distribution')
colorbar('vert')
axis equal

try
if results.sonicWarning==1
    figure(104)% 	Supersonic warning plot
    rotate3d on
    MAP=[0.7 1 0.7;1 0 0];
    colormap(MAP);
    h=fill3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)',results.sonicpanels');
    title('Supersonic Flow Warning. Red panels supersonic.')
    axis equal
end
end



%*************************************************
% Spanload plot, first attempt.

try
   figure(10)
   hold on, grid on
   plot(results.ystation(:,1),results.ForcePerMeter(:,1));
  title('Spanload on main wing');
  ylabel('Force per meter');
  xlabel('Spanstation')
catch
end
%Local CL plot
try
    
        figure(11)
        hold on, grid on
        plot(results.ystation(:,1),results.CL_local(:,1));
        title('Local CL on main wing');
        ylabel('CL');
        xlabel('Spanstation')
        text(0,0,strcat('Total CL ',num2str(results.CL)));
catch
end

try
   figure(12)
   hold on, grid on
   
   
   plot(results.ystation(:,1)./(ref.b_ref/2),results.CL_local(:,1)./results.CL);
   
   
   title('Normalized Local CL on main wing');
   
   ylabel('Normalized CL');
   xlabel('Normalized spanstation, \eta, [-]')
   
catch
end

try
   %% Bending moment diagram
   figure(13)
   hold on, grid on
   title('Bending moment on main wing');   
   ylabel('Bending moment, M_b, [Nm]');
   xlabel('Spanstation, y, [m]')
   plot(results.ystation,results.BendingMoment)
   


   figure(14)
    plot(results.ystation,results.ShearForce) ;
    hold on, grid on
   title('Shear force on main wing');   
   ylabel('Shear force, F_s, [N]');
   xlabel('Spanstation, y, [m]')
    
catch
    
  
end
%*************************************************
figure(7)
axis off
text(0,1,'Tornado Computation Results ')
text(0,.95,'JID: '); text(0.25,0.95,JID)
text(0.4,.95,'Downwash matrix condition: '); 
                                        text(0.8,0.95,num2str(results.dwcond))

text(0.4,.9,'Supersonic flow warning: ');
                                        try
                                        text(0.8,0.9,num2str(results.sonicWarning))
                                        catch
                                        text(0.8,0.9,'N/A')  
                                        end
                                        
text(0,.90,'Reference area: ');         text(0.25,0.90,num2str(ref.S_ref));	
text(0,.85,'Reference chord: ');        text(0.25,.85,num2str(ref.C_mac));
text(0,.8,'Reference span: ');          text(0.25,.8,num2str(ref.b_ref));  

text(0.4,.85,'Reference point pos: ');  text(0.8,.85,num2str(geo.ref_point));
try
text(0.4,.80,'Center of gravity  : ');  text(0.8,.80,num2str(geo.CG));
end
text(0,.7,'Net Wind Forces: (N)');			
   text(0.0,.65,'Drag: ');              text(0.1,.65,num2str(results.D));
   text(0.0,.6,'Side: ');           	text(0.1,.6,num2str(results.C));
   text(0.0,.55,'Lift: ');              text(0.1,.55,num2str(results.L));
   
text(0.35,.7,'Net Body Forces: (N)');			
   text(0.35,.65,'X: ');            	text(0.4,.65,num2str(results.FORCE(1)));
   text(0.35,.6,'Y: ');                 text(0.4,.6,num2str(results.FORCE(2)));
   text(0.35,.55,'Z: ');                text(0.4,.55,num2str(results.FORCE(3)));
    
text(0.7,.7,'Net Body Moments: (Nm)');
   text(0.7,.65,'Roll: ');              text(0.8,.65,num2str(results.MOMENTS(1)));
   text(0.7,.6,'Pitch: ');              text(0.8,.6,num2str(results.MOMENTS(2)));   
   text(0.7,.55,'Yaw: ');               text(0.8,.55,num2str(results.MOMENTS(3)));
   
text(0,.45,'CL ');                      text(0.1,.45,num2str(results.CL))   
text(0,.4, 'CD ');                      text(0.1,.4,num2str(results.CD)) 
text(0,.35,'CY ');                      text(0.1,.35,num2str(results.CY))
text(0,.30,'CD_t_r_e_f_f_t_z ');

try
    text(0.15,.30,num2str(results.Trefftz_drag_Coeff))
catch
    text(0.15,.30,'N/A')
end


text(0.35,.45,'CZ ');           text(0.45,.45,num2str(results.CZ))
text(0.35,.4,'CX ');            text(0.45,.4,num2str(results.CX))
text(0.35,.35,'CC ');           text(0.45,.35,num2str(results.CC))

text(0.7,.45,'Cm ');            text(0.8,.45,num2str(results.Cm))
text(0.7,.4,'Cn ');             text(0.8,.4,num2str(results.Cn))
text(0.7,.35,'Cl ');        	text(0.8,.35,num2str(results.Cl))

text(0,.2,'STATE: ');
text(0,.15,'\alpha [deg]: ');   text(.15,.15,num2str(state.alpha*180/pi));
text(0,.1,'\beta [deg]: ');     text(.15,.1,num2str(state.betha*180/pi));
text(0,.05,'Airspeed: ');       text(.15,.05,num2str(state.AS));
try
 text(0,.0,'Altitude: ');       text(.15,.0,num2str(state.ALT));
end
 text(0,-.05,'Density: ');      text(.15,-.05,num2str(state.rho));

text(0.3,.15,'P [rad/s]: ');    text(.45,.15,num2str(state.P));
text(0.3,.10,'Q [rad/s]: ');    text(.45,.10,num2str(state.Q));
text(0.3,.05,'R [rad/s]: ');    text(.45,.05,num2str(state.R));


    text(0.3,0.0,'PG Correction: ');
try
    text(0.55,0.0,num2str(state.pgcorr))
catch
    text(0.55,0,'N/A')
end
%text(0.6,.1,'Rudder setting [deg]:'); text(.9,.1,num2str(geo.flap_vector*180/pi));


[void sos voild]=ISAtmosphere(state.ALT);
Mach=state.AS/sos;
text(0.3,-.05,'Mach: ');    text(.45,-.05,num2str(Mach));

%%%%%
%
%%%%%
figure(8)
axis off
grid on
text(0,1,'TORNADO CALCULATION RESULTS, Derivatives')
text(0,.95,'JID: ');                text(0.25,0.95,JID)
text(0,.90,'Reference area: ');	    text(0.25,0.90,num2str( ref.S_ref ));	
text(0,.85,'Reference chord: ');    text(0.25,.85,num2str(ref.C_mac));
text(0,.8,'Reference span: ');      text(0.25,.8,num2str(ref.b_ref));
  
   
text(0.4,.90,'\alpha [deg]: '); 	text(.55,.9,num2str(state.alpha*180/pi));
text(0.4,.85, '\beta [deg]: ');     text(.55,.85,num2str(state.betha*180/pi));
text(0.4,.8,'Airspeed: ');          text(.55,.8,num2str(state.AS));

text(0.65,.9, 'P [rad/s]: ');  		text(.8,.9,num2str(state.P));
text(0.65,.85,'Q[rad/s]: '); 		text(.8,.85,num2str(state.Q));
text(0.65,.8,' R[rad/s]: ');  		text(.8,.8,num2str(state.R));

figure(8)
text(0,.7,'CL derivatives : ');			

text(0,.65,'CL_{\alpha}');		    text(0.15,.65,num2str(results.CL_a));
text(0,.6,'CL_{\beta}');		    text(0.15,.6,num2str(results.CL_b));
text(0,.55,'CL_P');			        text(0.15,.55,num2str(results.CL_P));
text(0,.5,'CL_Q');			        text(0.15,.5,num2str(results.CL_Q));
text(0,.45,'CL_R');			        text(0.15,.45,num2str(results.CL_R));
   
text(0,.35,'Roll derivatives : ');			
   
text(0,.3,'Cl_{\alpha}');		    text(0.15,.3,num2str(results.Cl_a));
text(0,.25,'Cl_{\beta}');		    text(0.15,.25,num2str(results.Cl_b));
text(0,.2,'Cl_P');			        text(0.15,.2,num2str(results.Cl_P));
text(0,.15,'Cl_Q');			        text(0.15,.15,num2str(results.Cl_Q));
text(0,.1,'Cl_R');			        text(0.15,.1,num2str(results.Cl_R));
   
text(0.35,.7,'CD derivatives : ');			

text(0.35,.65,'CD_{\alpha}');	text(0.5,.65,num2str(results.CD_a));
text(0.35,.6,'CD_{\beta}');	text(0.5,.6,num2str(results.CD_b));
text(0.35,.55,'CD_P');		text(0.5,.55,num2str(results.CD_P));
text(0.35,.5,'CD_Q');		text(0.5,.5,num2str(results.CD_Q));
text(0.35,.45,'CD_R');		text(0.5,.45,num2str(results.CD_R));
   
text(0.35,.35,'Pitch derivatives : ');			

text(0.35,.3,'Cm_{\alpha}');	text(.5,.3,num2str(results.Cm_a));
text(0.35,.25,'Cm_{\beta}');	text(0.5,.25,num2str(results.Cm_b));
text(0.35,.2,'Cm_P');		text(0.5,.2,num2str(results.Cm_P));
text(0.35,.15,'Cm_Q');		text(0.5,.15,num2str(results.Cm_Q));
text(0.35,.1,'Cm_R');		text(0.5,.1,num2str(results.Cm_R));
   
text(0.7,.7,'CY derivatives : ');			

text(0.7,.65,'CY_{\alpha}');	text(0.85,.65,num2str(results.CY_a));
text(0.7,.6,'CY_{\beta}');		text(0.85,.6,num2str(results.CY_b));
text(0.7,.55,'CY_P');		text(0.85,.55,num2str(results.CY_P));
text(0.7,.5,'CY_Q');			text(0.85,.5,num2str(results.CY_Q));
text(0.7,.45,'CY_R');		text(0.85,.45,num2str(results.CY_R));
   
text(0.7,.35,'Yaw derivatives : ');			

text(0.7,.3,'Cn_{\alpha}');	text(.85,.3,num2str(results.Cn_a));
text(0.7,.25,'Cn_{\beta}');	text(0.85,.25,num2str(results.Cn_b));
text(0.7,.2,'Cn_P');			text(0.85,.2,num2str(results.Cn_P));
text(0.7,.15,'Cn_Q');		text(0.85,.15,num2str(results.Cn_Q));
text(0.7,.1,'Cn_R');			text(0.85,.1,num2str(results.Cn_R));
   
   
figure(9)

axis off
%grid on
text(0,1,'TORNADO CALCULATION RESULTS, Central difference, RUDDER DERIVs')
text(0,.95,'JID: '); text(0.25,0.95,JID)

text(0,.90,'Reference area: ');	text(0.25,0.90,num2str(ref.S_ref));	
text(0,.85,'Reference chord: ');text(0.25,.85,num2str(ref.C_mac));
text(0,.8,'Reference span: ');text(0.25,.8,num2str(ref.b_ref));
   
text(0.4,.90,'\alpha: '); 	text(.55,.9,num2str(state.alpha*180/pi));
text(0.4,.85,'\beta: ');     text(.55,.85,num2str(state.betha*180/pi));
text(0.4,.8,'Airspeed: ');  text(.55,.8,num2str(state.AS));

text(0.65,.9,'P: ');  		text(.8,.9,num2str(state.P));
text(0.65,.85,'Q: '); 		text(.8,.85,num2str(state.Q));
text(0.65,.8,'R: ');  		text(.8,.8,num2str(state.R));

text(0,.45,'CL_{\delta}');			text(0.15,.45,num2str(results.CL_d'));
text(0,.05,'Cl_{\delta}');			text(0.15,.05,num2str(results.Cl_d'));
text(0.35,.45,'CD_{\delta}');		text(0.5,.45,num2str(results.CD_d'));
text(0.35,.05,'Cm_{\delta}');		text(0.5,.05,num2str(results.Cm_d'));
text(0.7,.45,'CY_{\delta}');		text(0.85,.45,num2str(results.CY_d'));
text(0.7,.05,'Cn_{\delta}');		text(0.85,.05,num2str(results.Cn_d'));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,z]=midpoint(XYZ);

[a b c]=size(XYZ);

for i=1:a
   x(i)=((XYZ(i,1,1)+XYZ(i,2,1))*2+XYZ(i,3,1)+XYZ(i,4,1))/6;
   y(i)=((XYZ(i,1,2)+XYZ(i,2,2))+XYZ(i,3,2)+XYZ(i,4,2))/4;
   z(i)=((XYZ(i,1,3)+XYZ(i,2,3))+XYZ(i,3,3)+XYZ(i,4,3))/4;
end
end%function
