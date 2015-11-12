function [bc,bci]=setDynamicBoundary(lattice,state,geo,h)
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
% usage: [boundarycondition] = boundary4(lattice,state,geo)
%
% This function computes the right hand side of the vortex lattice equation
% system. I.e. the velocity parallell to the panel normal through each 
% collocation point due to rotations and angle of attack and sideslip.
%
% Example:
%
%  rhs=(setboundary4(lattice,state,geo))';
%
% Calls:
%           None    
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado core function
%
% Revision History:
%   Cramer,   2015 06 23:  Adapted for array alpha and airspeed
%   Cramer,   2015 05 07:  Removed flapped parts
%   Bristol,  2007 06 27:  Addition of new header. TM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b c]=size(lattice.COLLOC);

try
    geo.CG;
catch
    geo.CG=geo.ref_point;
end


%%%%
%Steady state boundary condition column
V=state.AS;

Wind=([V.*cos(state.alpha).*cos(state.betha) -V.*cos(state.alpha).*sin(state.betha) V.*sin(state.alpha)]);

for i=1:a
   Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]);
end     

veloc=Wind+Rot;
bc=sum(lattice.N.*veloc,2)';  %steady state bc
%%%%%%%%%%%%%%%%%%%%%%%

%%%%
%Imaginary boundary condition airspeed derivative
V=state.AS+1j*h;

Wind=([V.*cos(state.alpha).*cos(state.betha) -V.*cos(state.alpha).*sin(state.betha) V.*sin(state.alpha)]);

for i=1:a
   Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]);
end     

veloc=Wind+Rot;
bci(:,1)=sum(lattice.N.*veloc,2)';  %airspeed derivative
%%%%%%%%%%%%%%%%%%%%%%%

%%%%
%Imaginary boundary condition angle of attack derivative
V=state.AS;

Wind=([V.*cos(state.alpha+1j*h).*cos(state.betha) -V.*cos(state.alpha+1j*h).*sin(state.betha) V.*sin(state.alpha+1j*h)]);

for i=1:a
   Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]);
end     

veloc=Wind+Rot;
bci(:,2)=sum(lattice.N.*veloc,2)';  %angle of attack derivative
%%%%%%%%%%%%%%%%%%%%%%%

%%%%
%Imaginary boundary condition beta derivative
V=state.AS;

Wind=([V.*cos(state.alpha).*cos(state.betha+1j*h) -V.*cos(state.alpha).*sin(state.betha+1j*h) V.*sin(state.alpha)]);

for i=1:a
   Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]);
end     

veloc=Wind+Rot;
bci(:,3)=sum(lattice.N.*veloc,2)';  %beta derivative
%%%%%%%%%%%%%%%%%%%%%%%

end %function