function Error = findCG2015_8_31(CG,Results)

global M C K SegNum halfWingLength c ActLoc cordNum bw numSpanB rho_inf

[masterStateQ2,masterLatticeQ2,masterGeoQ2,masterResultsQ2,refQ2,forceTwistQ2,CLQ2,CDQ2] = staticAsymmetricAlphaSweep2015_8_29(M,C,K,SegNum,halfWingLength,c,CG,ActLoc,cordNum,bw,c,0,0,numSpanB,[-4:10],[-6:6],0,2,rho_inf);
[masterStateQ3,masterLatticeQ3,masterGeoQ3,masterResultsQ3,refQ3,forceTwistQ3,CLQ3,CDQ3] = staticAsymmetricAlphaSweep2015_8_29(M,C,K,SegNum,halfWingLength,c,CG,ActLoc,cordNum,bw,c,0,0,numSpanB,[-4:10],[-6:6],0,3,rho_inf);
[masterStateQ4,masterLatticeQ4,masterGeoQ4,masterResultsQ4,refQ4,forceTwistQ4,CLQ4,CDQ4] = staticAsymmetricAlphaSweep2015_8_29(M,C,K,SegNum,halfWingLength,c,CG,ActLoc,cordNum,bw,c,0,0,numSpanB,[-4:10],[-6:6],0,4,rho_inf);

SimResults = [[masterResultsQ2(1,1:2:end).Cl]';[masterResultsQ2(1,1:2:end).Cm]';[masterResultsQ2(1,1:2:end).Cn]';[masterResultsQ2(5,1:2:end).Cl]';[masterResultsQ2(5,1:2:end).Cm]';[masterResultsQ2(5,1:2:end).Cn]';...
    [masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';[masterResultsQ2(7,1:2:end).Cl]';[masterResultsQ2(7,1:2:end).Cm]';[masterResultsQ2(7,1:2:end).Cn]';...
    [masterResultsQ2(9,1:2:end).Cl]';[masterResultsQ2(9,1:2:end).Cm]';[masterResultsQ2(9,1:2:end).Cn]';[masterResultsQ2(11,1:2:end).Cl]';[masterResultsQ2(11,1:2:end).Cm]';[masterResultsQ2(11,1:2:end).Cn]';[masterResultsQ2(13,1:2:end).Cl]';[masterResultsQ2(13,1:2:end).Cm]';[masterResultsQ2(13,1:2:end).Cn]';...
    [masterResultsQ3(7,1:2:end).Cl]';[masterResultsQ3(7,1:2:end).Cm]';[masterResultsQ3(7,1:2:end).Cn]';[masterResultsQ3(9,1:2:end).Cl]';[masterResultsQ3(9,1:2:end).Cm]';[masterResultsQ3(9,1:2:end).Cn]';[masterResultsQ3(11,1:2:end).Cl]';[masterResultsQ3(11,1:2:end).Cm]';[masterResultsQ3(11,1:2:end).Cn]';[masterResultsQ3(13,1:2:end).Cl]';[masterResultsQ3(13,1:2:end).Cm]';[masterResultsQ3(13,1:2:end).Cn]';...
    [masterResultsQ4(7,1:2:end).Cl]';[masterResultsQ4(7,1:2:end).Cm]';[masterResultsQ4(7,1:2:end).Cn]';[masterResultsQ4(9,1:2:end).Cl]';[masterResultsQ4(9,1:2:end).Cm]';[masterResultsQ4(9,1:2:end).Cn]';[masterResultsQ4(11,1:2:end).Cl]';[masterResultsQ4(11,1:2:end).Cm]';[masterResultsQ4(11,1:2:end).Cn]';[masterResultsQ4(13,1:2:end).Cl]';[masterResultsQ4(13,1:2:end).Cm]';[masterResultsQ4(13,1:2:end).Cn]'];
Error = Results-SimResults;
temp = [];
for i = 1:length(Error)/24
    temp = [temp;Error((i-1)*24+1:(i-1)*24+9);Error((i-1)*24+17:(i-1)*24+24)];
end
Error = norm(temp);

disp(['CG: ',num2str(CG(1)),', ',num2str(CG(2))])
disp(['Error: ',num2str(Error)])