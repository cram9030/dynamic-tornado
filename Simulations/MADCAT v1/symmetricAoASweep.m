function [spanCL,Lift,Drag,Pitching,Results,Lattice,State] = symmetricAoASweep(geo,S_ref,C_ref,B_ref,mac_pos,AoA,speed,airDensity,totalPanels,twist)

%Intialize output
spanCL = zeros(geo.ny,length(AoA));
Pitching = zeros(length(AoA),1);
Lift = zeros(length(AoA),1);
Drag = zeros(length(AoA),1);

%Set refrence values for coefficients
ref=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos);
%Set the perscribed twist
if twist == 0
    geo.Wings.wing.Theta = zeros(size(geo.Wings.wing.Theta));
else
    geo.Wings.wing.Theta = [flipud([0:twist/(geo.ny/2):twist]');[twist/(geo.ny/2):twist/(geo.ny/2):twist]']*pi/180;
end

for i = 1:length(AoA)
    %Create state structure reflecting the new angle of attack
    state = setupState2015_8_10(AoA(i)*pi/180*ones(totalPanels,1),AoA(i)*pi/180,0,speed*ones(totalPanels,1),speed,airDensity);
    %Generate the lattice for VLM calulations
    lattice = generateLattice(geo,S_ref,C_ref,B_ref,mac_pos,state);
    %Solve the VLM
    results=dynamicSolver(state,geo,lattice,lattice,1e-20,0);
    %Calculate the coefficients using the refrence values
    results = coeff_create(results,lattice,state,ref,geo,0);
    spanCL(:,i) = results.spanResults.CL_local;
    Pitching(i) = results.MOMENTS(2);
    Lift(i) = results.L;
    Drag(i) = results.D;
    Results(i) = results;
    Lattice(i) = lattice;
    State(i) = state;
end