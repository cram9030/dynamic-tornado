function [Lift,Drag,Side,Pitching,Roll,Yaw] = asymmetricAoASweep(geo,S_ref,C_ref,B_ref,mac_pos,AoA,speed,airDensity,totalPanels)

%Initalize output parameters
Pitching = zeros(length(AoA),1);
Roll = zeros(length(AoA),1);
Yaw = zeros(length(AoA),1);
Lift = zeros(length(AoA),1);
Drag = zeros(length(AoA),1);
Side = zeros(length(AoA),1);

%Set refrence values for coefficients
ref=setRef2015_8_17(S_ref,C_ref,B_ref,mac_pos);
for i = 1:length(AoA)
    %Create state structure reflecting the new angle of attack
    state = setupState2015_8_10(AoA(i)*pi/180*ones(totalPanels,1),AoA(i)*pi/180,0,speed*ones(totalPanels,1),speed,airDensity);
    %Generate the lattice for VLM calulations
    lattice = generateLattice(geo,S_ref,C_ref,B_ref,mac_pos,state);
    %Solve the VLM
    results=dynamicSolver(state,geo,lattice,lattice,1e-20,0);
    %Calculate the coefficients using the refrence values
    results = coeff_create(results,lattice,state,ref,geo,0);
    %Populate moments
    Pitching(i) = results.MOMENTS(2);
    Roll(i) = results.MOMENTS(1);
    Yaw(i) = results.MOMENTS(3);
    Lift(i) = results.L;
    Drag(i) = results.D;
    Side(i) = results.C;
end