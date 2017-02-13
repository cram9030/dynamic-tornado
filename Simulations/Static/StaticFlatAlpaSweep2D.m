clear masterState
clear masterResults

CLNoBody = zeros(15,1);
CDNoBody = zeros(15,1);

[ref]=setRef(L,bw,c,c,centroid,0);

for i = -4:10
    masterState(i+5) = setupState5_7_2015(i,0,U_inf,rhom_inf);
    [lattice] = generateLattice5_6(SegNum,zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),c,centroid,cordNum,L,bw,c,0,0,numSpanB,masterState(i+5));
    [geo] = setupGeo5_7_2015(zeros(1,3),zeros(1,3),2*max(max(max(lattice.XYZ))),cordNum,SegNum,numSpanB);
    results = solver5_7_2015(masterState(i+5),geo,lattice);
    masterResults(i+5) = coeff_create5_7_2015(results,lattice,masterState(i+5),ref,geo);
    CLNoBody(i+5) = masterResults(i+5).CL;
    CDNoBody(i+5) = masterResults(i+5).CD;
    disp(['Precent Complete: ',num2str((i+5)/length(-4:10))])
end

plot(-4:10,CLNoBody)
figure;plot(-4:10,CDNoBody)

% clear masterState
% clear masterResults
% 
% CLNoBiPlane = zeros(15,1);
% CDNoBiPlane = zeros(15,1);
% 
% [ref]=setRef(L,bw,c,c,centroid,0);
% 
% for i = -4:10
%     masterState(i+5) = setupState5_7_2015(i,0,U_inf,rhom_inf);
%     [lattice] = generateLattice5_27(SegNum,zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),c,centroid,cordNum,L,bw,bc,noseH,zb,numSpanB,masterState(i+5));
%     %[lattice] = generateLattice5_18(SegNum,zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),zeros(SegNum+1,1),c,centroid,cordNum,L,bw,c,0,0,numSpanB,masterState(i+5));
%     [geo] = setupGeo5_7_2015(zeros(1,3),zeros(1,3),2*max(max(max(lattice.XYZ))),cordNum,SegNum,numSpanB);
%     results = solver5_7_2015(masterState(i+5),geo,lattice);
%     masterResults(i+5) = coeff_create5_7_2015(results,lattice,masterState(i+5),ref,geo);
%     CLNoBiPlane(i+5) = masterResults(i+5).CL;
%     CDNoBiPlane(i+5) = masterResults(i+5).CD;
%     disp(['Precent Complete: ',num2str((i+5)/length(-4:10))])
% end
% 
% plot(-4:10,CLNoBiPlane)
% figure;plot(-4:10,CDNoBiPlane)