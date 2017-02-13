clear masterState
clear masterResults

CL = zeros(15,1);
CD = zeros(15,1);

[ref]=setRef(L,bw,c,c,centroid,0);

for i = -4:10
    disp(['Precent Complete: ',num2str((i+5)/length(-4:10))])
    masterState(i+5) = setupState5_7_2015(i,0,U_inf,rhom_inf);
    [lattice] = generateLattice5_21(SegNum,zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),c,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,masterState(i+5));
    [geo] = setupGeo5_7_2015(zeros(1,3),zeros(1,3),2*max(max(max(lattice.XYZ))),cordNum,SegNum,numSpanB);
    results = solver5_7_2015(masterState(i+5),geo,lattice);
    masterResults(i+5) = coeff_create5_7_2015(results,lattice,masterState(i+5),ref,geo);
    CL(i+5) = masterResults(i+5).CL;
    CD(i+5) = masterResults(i+5).CD;
end

plot(-4:10,CL)
figure;plot(-4:10,CD)