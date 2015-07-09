function [xNose,yNose,zNose,xNoseLattice,yNoseLattice,zNoseLattice,xNosem,yNosem,zNosem,sNose] = generateNose5_27(numSpanB,bw,centroid,noseH,zb,cordNum,cord)

xNose = [];
yNose = [];
zNose = [];
xNoseLattice = [];
yNoseLattice = [];
zNoseLattice = [];
xNosem = [];
yNosem = [];
zNosem = [];
sNose = [];

numMid = ceil((0.0254/4)/(0.038989/numSpanB/2));
if mod(numMid,2) == 1
    numMid = numMid + 1;
end

yi = [-bw/2:bw/numSpanB:bw/2];

xi = sqrt(noseH^2-yi.^2*noseH^2/(bw/2)^2)+centroid(1)+noseH;
xi(1) = centroid(1)+noseH;
xi(end) = centroid(1)+noseH;

zi = sqrt(0.038989^2/4-0.038989^2/4*(xi-centroid(1)-noseH).^2/noseH^2) + 0.0254/4;
zi(ceil(length(yi)/2)) = 0.0254/4;

xi = xi';
yi = yi';
zi = zi';

xNose = [[xi(1)-(cord-noseH)/cordNum;xi(1);xi(2);xi(1)-(cord-noseH)/cordNum;xi(1)-(cord-noseH)/cordNum],...
    [xi(end)-(cord-noseH)/cordNum;xi(end);xi(end-1);xi(end)-(cord-noseH)/cordNum;xi(end)-(cord-noseH)/cordNum]];
yNose = [[yi(1);yi(1);yi(2);yi(2);yi(1)],...
    [yi(end);yi(end);yi(end-1);yi(end-1);yi(end)]];
zNose = [[zi(1);zi(1);zi(2);zi(1);zi(1)],...
    [zi(end);zi(end);zi(end-1);zi(end);zi(end)]];
xNoseLattice = [[xi(1)-(cord-noseH)/cordNum;xi(1)-.25*(cord-noseH)/cordNum;(-xi(2)+(xi(1)-(cord-noseH)/cordNum))/4+xi(2);xi(1)-(cord-noseH)/cordNum],...
    [xi(end)-(cord-noseH)/cordNum;xi(end)-.25*(cord-noseH)/cordNum;(-xi(end-1)+(xi(end)-(cord-noseH)/cordNum))/4+xi(end-1);xi(end)-(cord-noseH)/cordNum]];
yNoseLattice = [[yi(1);yi(1);yi(2);yi(2)],...
    [yi(end);yi(end);yi(end-1);yi(end-1)]];
zNoseLattice = [[zi(1);zi(1);(zi(1)-zi(2))/4+zi(2);zi(1)],...
    [zi(end);zi(end);(zi(end)-zi(end-1))/4+zi(end-1);zi(end)]];
xNosem = [((xi(2)+xi(1))/2-xi(1)+(cord-noseH)/cordNum)/4+xi(1)-(cord-noseH)/cordNum,((xi(end)+xi(end-1))/2-xi(end)+(cord-noseH)/cordNum)/4+xi(end)-(cord-noseH)/cordNum];
yNosem = [(yi(1)+yi(2))/2,(yi(end)+yi(end-1))/2];
zNosem = [(-zi(1)+zi(2))/8+zi(1),(-zi(end)+zi(end-1))/8+zi(end)];
sNose = [pi,0];

for i = 2:numSpanB/2
    xNose = [xNose,[xi(1)-(cord-noseH)/cordNum;xi(1);xi(1);xi(1)-(cord-noseH)/cordNum;xi(1)-(cord-noseH)/cordNum],...
        [xi(end)-(cord-noseH)/cordNum;xi(end);xi(end);xi(end)-(cord-noseH)/cordNum;xi(end)-(cord-noseH)/cordNum]];
    yNose = [yNose,[yi(i);yi(i);yi(i+1);yi(i+1);yi(i)],...
        [yi(end-i+1);yi(end-i+1);yi(end-i);yi(end-i);yi(end-i+1)]];
    zNose = [zNose,[zi(1);zi(1);zi(1);zi(1);zi(1)],...
        [zi(1);zi(1);zi(1);zi(1);zi(1)]];
    xNoseLattice = [xNoseLattice,[xi(1)-(cord-noseH)/cordNum;xi(1)-.25*(cord-noseH)/cordNum;xi(1)-.25*(cord-noseH)/cordNum;xi(1)-(cord-noseH)/cordNum],...
        [xi(end)-(cord-noseH)/cordNum;xi(end)-.25*(cord-noseH)/cordNum;xi(end)-.25*(cord-noseH)/cordNum;xi(end)-(cord-noseH)/cordNum]];
    yNoseLattice = [yNoseLattice,[yi(i);yi(i);yi(i+1);yi(i+1)],...
        [yi(end-i+1);yi(end-i+1);yi(end-i);yi(end-i)]];
    zNoseLattice = [zNoseLattice,[zi(1);zi(1);zi(1);zi(1)],...
        [zi(1);zi(1);zi(1);zi(1)]];
    xNosem = [xNosem,xi(1)-.75*(cord-noseH)/cordNum,xi(end)-.75*(cord-noseH)/cordNum];
    yNosem = [yNosem,(yi(i)+yi(i+1))/2,(yi(end-i+1)+yi(end-i))/2];
    zNosem = [zNosem,zi(1),zi(1)];
    sNose = [sNose,0,pi];
end

%Create top left edges
for i = 3:ceil(length(yi)/2)
    c = xi(i-1)-xi(i-2)+(xi(i)-xi(i-1))/2;
    xNose = [xNose,[xi(i-2);xi(i-1);xi(i);xi(i-2);xi(i-2)]];
    yNose = [yNose,[yi(i-1);yi(i-1);yi(i);yi(i);yi(i-1)]];
    zNose = [zNose,[zi(i-2);zi(i-1);zi(i);zi(i-2);zi(i-2)]];
    xNoseLattice = [xNoseLattice,[xi(i-2);xi(i-1)-c/4;xi(i)-c/4;xi(i-2)]];
    yNoseLattice = [yNoseLattice,[yi(i-1);yi(i-1);yi(i);yi(i)]];
    zNoseLattice = [zNoseLattice,[zi(i-2);zi(i-1)+(zi(i-2)-zi(i-1))/4;zi(i)+(zi(i-2)-zi(i))/4;zi(i-2)]];
    xNosem = [xNosem,xi(i-2)+c/4];
    yNosem = [yNosem,(yi(i)-yi(i-1))/2+yi(i-1)];
    zNosem = [zNosem,zi(i-2)-(zi(i-2)-(zi(i-1)-zi(i))/2-zi(i))/4];
    sNose = [sNose,pi];
end

%Create top right edges
for i = ceil(length(yi)/2):length(yi)-2
    c = xi(i)-xi(i+2)+(xi(i+1)-xi(i))/2;
    xNose = [xNose,[xi(i+2);xi(i+1);xi(i);xi(i+2);xi(i+2)]];
    yNose = [yNose,[yi(i+1);yi(i+1);yi(i);yi(i);yi(i+1)]];
    zNose = [zNose,[zi(i+2);zi(i+1);zi(i);zi(i+2);zi(i+2)]];
    xNoseLattice = [xNoseLattice,[xi(i+2);xi(i+1)-c/4;xi(i)-c/4;xi(i+2)]];
    yNoseLattice = [yNoseLattice,[yi(i+1);yi(i+1);yi(i);yi(i)]];
    zNoseLattice = [zNoseLattice,[zi(i+2);zi(i+1)+(zi(i+2)-zi(i+1))/4;zi(i)+(zi(i+2)-zi(i))/4;zi(i+2)]];
    xNosem = [xNosem,xi(i+2)+c/4];
    yNosem = [yNosem,(yi(i+1)-yi(i))/2+yi(i)];
    zNosem = [zNosem,zi(i+2)-(zi(i+2)-(zi(i+1)-zi(i))/2-zi(i))/4];
    sNose = [sNose,0];
end

% %Create top left edges
% for i = 3:ceil(length(yi)/2)
%     c = xi(i-1)-xi(i-2)+(xi(i)-xi(i-1))/2;
%     xNose = [xNose,[xi(i-2);xi(i-1);xi(i);xi(i-2);xi(i-2)],[xi(i-2);xi(i-1);xi(i);xi(i-2);xi(i-2)]];
%     yNose = [yNose,[yi(i-1);yi(i-1);yi(i);yi(i);yi(i-1)],[yi(i-2);yi(i-1);yi(i);yi(i-2);yi(i-2)]];
%     zNose = [zNose,[zi(i-2);zi(i-1);zi(i);zi(i-2);zi(i-2)],[zi(i-1);zi(i-1);zi(i);zi(i);zi(i-1)]];
%     xNoseLattice = [xNoseLattice,[xi(i-2);xi(i-1)-c/4;xi(i)-c/4;xi(i-2)],[xi(i-2);xi(i-1)-c/4;xi(i)-c/4;xi(i-2)]];
%     yNoseLattice = [yNoseLattice,[yi(i-1);yi(i-1);yi(i);yi(i)],[yi(i-2);yi(i-1)+(yi(i-2)-yi(i-1))/4;yi(i)+(yi(i-2)-yi(i))/4;yi(i-2)]];
%     zNoseLattice = [zNoseLattice,[zi(i-2);zi(i-1)+(zi(i-2)-zi(i-1))/4;zi(i)+(zi(i-2)-zi(i))/4;zi(i-2)],[zi(i-1);zi(i-1);zi(i);zi(i)]];
%     xNosem = [xNosem,xi(i-2)+c/4,xi(i-2)+c/4];
%     yNosem = [yNosem,(yi(i)-yi(i-1))/2+yi(i-1),yi(i-2)-(yi(i-2)-(yi(i-1)-yi(i))/2-yi(i))/4];
%     zNosem = [zNosem,zi(i-2)-(zi(i-2)-(zi(i-1)-zi(i))/2-zi(i))/4,(zi(i)-zi(i-1))/2+zi(i-1)];
%     sNose = [sNose,0,pi];
% end
% 
% %Create top right edges
% for i = ceil(length(yi)/2):length(yi)-2
%     c = xi(i)-xi(i+2)+(xi(i+1)-xi(i))/2;
%     xNose = [xNose,[xi(i+2);xi(i+1);xi(i);xi(i+2);xi(i+2)],[xi(i+2);xi(i+1);xi(i);xi(i+2);xi(i+2)]];
%     yNose = [yNose,[yi(i+1);yi(i+1);yi(i);yi(i);yi(i+1)],[yi(i+2);yi(i+1);yi(i);yi(i+2);yi(i+2)]];
%     zNose = [zNose,[zi(i+2);zi(i+1);zi(i);zi(i+2);zi(i+2)],[zi(i+1);zi(i+1);zi(i);zi(i);zi(i+1)]];
%     xNoseLattice = [xNoseLattice,[xi(i+2);xi(i+1)-c/4;xi(i)-c/4;xi(i+2)],[xi(i+2);xi(i+1)-c/4;xi(i)-c/4;xi(i+2)]];
%     yNoseLattice = [yNoseLattice,[yi(i+1);yi(i+1);yi(i);yi(i)],[yi(i+2);yi(i+1)+(yi(i+2)-yi(i+1))/4;yi(i)+(yi(i+2)-yi(i))/4;yi(i+2)]];
%     zNoseLattice = [zNoseLattice,[zi(i+2);zi(i+1)+(zi(i+2)-zi(i+1))/4;zi(i)+(zi(i+2)-zi(i))/4;zi(i+2)],[zi(i+1);zi(i+1);zi(i);zi(i)]];
%     xNosem = [xNosem,xi(i+2)+c/4,xi(i+2)+c/4];
%     yNosem = [yNosem,(yi(i+1)-yi(i))/2+yi(i),yi(i+2)-(yi(i+2)-(yi(i+1)-yi(i))/2-yi(i))/4];
%     zNosem = [zNosem,zi(i+2)-(zi(i+2)-(zi(i+1)-zi(i))/2-zi(i))/4,(zi(i+1)-zi(i))/2+zi(i)];
%     sNose = [sNose,0,pi];
% end

%Left Top
for i = 1:ceil(length(yi)/2)-2
    for j = i+2:ceil(length(yi)/2)-1
        c = xi(i+1)-xi(i);
        xNose = [xNose,[xi(i);xi(i+1);xi(i+1);xi(i);xi(i)]];
        yNose = [yNose,[yi(j);yi(j);yi(j+1);yi(j+1);yi(j)]];
        zNose = [zNose,[zi(i);zi(i+1);zi(i+1);zi(i);zi(i)]];
        xNoseLattice = [xNoseLattice,[xi(i);xi(i+1)-c/4;xi(i+1)-c/4;xi(i)]];
        yNoseLattice = [yNoseLattice,[yi(j);yi(j);yi(j+1);yi(j+1)]];
        zNoseLattice = [zNoseLattice,[zi(i);zi(i+1)-(zi(i+1)-zi(i))/4;zi(i+1)-(zi(i+1)-zi(i))/4;zi(i)]];
        xNosem = [xNosem,xi(i)+c/4];
        yNosem = [yNosem,(yi(j+1)-yi(j))/2+yi(j)];
        zNosem = [zNosem,(zi(i+1)-zi(i))/4+zi(i)];
        sNose = [sNose,0];
    end
end

% %Create Left side
% for i = 1:ceil(length(yi)/2)-3
%     c = xi(i+1)-xi(i);
%     for j = i+2:ceil(length(yi)/2)-1
%         xNose = [xNose,[xi(i);xi(i+1);xi(i+1);xi(i);xi(i)]];
%         yNose = [yNose,[yi(i);yi(i+1);yi(i+1);yi(i);yi(i)]];
%         zNose = [zNose,[zi(j);zi(j);zi(j+1);zi(j+1);zi(j)]];
%         xNoseLattice = [xNoseLattice,[xi(i);xi(i+1)-c/4;xi(i+1)-c/4;xi(i)]];
%         yNoseLattice = [yNoseLattice,[yi(i);yi(i+1)-(yi(i+1)-yi(i))/4;yi(i+1)-(yi(i+1)-yi(i))/4;yi(i)]];
%         zNoseLattice = [zNoseLattice,[zi(j);zi(j);zi(j+1);zi(j+1)]];
%         xNosem = [xNosem,xi(i)+c/4];
%         yNosem = [yNosem,(yi(i+1)-yi(i))/4+yi(i)];
%         zNosem = [zNosem,(zi(j+1)-zi(j))/2+zi(j)];
%         sNose = [sNose,pi];
%     end
% end
% 
% for i = 1:ceil(length(yi)/2)-1
%     c = xi(i+1)-xi(i);
%     for j = 0:numMid/2-1
%         xNose = [xNose,[xi(i);xi(i+1);xi(i+1);xi(i);xi(i)]];
%         yNose = [yNose,[yi(i);yi(i+1);yi(i+1);yi(i);yi(i)]];
%         zNose = [zNose,[j*0.0254/numMid;j*0.0254/numMid;(j+1)*0.0254/numMid;(j+1)*0.0254/numMid;j*0.0254/numMid]];
%         xNoseLattice = [xNoseLattice,[xi(i);xi(i+1)-c/4;xi(i+1)-c/4;xi(i)]];
%         yNoseLattice = [yNoseLattice,[yi(i);yi(i+1)-(yi(i+1)-yi(i))/4;yi(i+1)-(yi(i+1)-yi(i))/4;yi(i)]];
%         zNoseLattice = [zNoseLattice,[j*0.0254/numMid;j*0.0254/numMid;(j+1)*0.0254/numMid;(j+1)*0.0254/numMid]];
%         xNosem = [xNosem,xi(i)+c/4];
%         yNosem = [yNosem,(yi(i+1)-yi(i))/4+yi(i)];
%         zNosem = [zNosem,(j+.5)*0.0254/numMid];
%         sNose = [sNose,0];
%     end
% end

%Create Right Top
for i = ceil(length(yi)/2)+1:length(yi)-2
    for j = i:length(yi)-2
        c = xi(j+2)-xi(j+1);
        xNose = [xNose,[xi(j+1);xi(j+2);xi(j+2);xi(j+1);xi(j+1)]];
        yNose = [yNose,[yi(i);yi(i);yi(i-1);yi(i-1);yi(i)]];
        zNose = [zNose,[zi(j+1);zi(j+2);zi(j+2);zi(j+1);zi(j+1)]];
        xNoseLattice = [xNoseLattice,[xi(j+1);xi(j+2)-3*c/4;xi(j+2)-3*c/4;xi(j+1)]];
        yNoseLattice = [yNoseLattice,[yi(i);yi(i);yi(i-1);yi(i-1)]];
        zNoseLattice = [zNoseLattice,[zi(j+1);zi(j+2)-3*(zi(j+2)-zi(j+1))/4;zi(j+2)-3*(zi(j+2)-zi(j+1))/4;zi(j+1)]];
        xNosem = [xNosem,xi(j+1)+3*c/4];
        yNosem = [yNosem,yi(i-1)+(yi(i)-yi(i-1))/2];
        zNosem = [zNosem,zi(j+1)+(zi(j+2)-zi(j+1))/2];
        sNose = [sNose,pi];
    end
end

% %Create Right Side
% for i = 0:ceil(length(yi)/2)-3
%     c = xi(end-i-1)-xi(end-i);
%     for j = i+2:ceil(length(yi)/2)-1
%         xNose = [xNose,[xi(end-i);xi(end-i-1);xi(end-i-1);xi(end-i);xi(end-i)]];
%         yNose = [yNose,[yi(end-i);yi(end-i-1);yi(end-i-1);yi(end-i);yi(end-i)]];
%         zNose = [zNose,[zi(j);zi(j);zi(j+1);zi(j+1);zi(j)]];
%         xNoseLattice = [xNoseLattice,[xi(end-i);xi(end-i-1);xi(end-i-1);xi(end-i)]];
%         yNoseLattice = [yNoseLattice,[yi(end-i);yi(end-i-1)-(yi(end-i-1)-yi(end-i))/4;yi(end-i-1)-(yi(end-i-1)-yi(end-i))/4;yi(end-i)]];
%         zNoseLattice = [zNoseLattice,[zi(j);zi(j);zi(j+1);zi(j+1)]];
%         xNosem = [xNosem,xi(end-i)+c/4];
%         yNosem = [yNosem,(yi(end-i-1)-yi(end-i))/4+yi(end-i)];
%         zNosem = [zNosem,(zi(j+1)-zi(j))/2+zi(j)];
%         sNose = [sNose,0];
%     end
% end
% 
% for i = 0:ceil(length(yi)/2)-2
%     c = xi(end-i-1)-xi(end-i);
%     for j = 0:numMid/2-1
%         xNose = [xNose,[xi(end-i);xi(end-i-1);xi(end-i-1);xi(end-i);xi(end-i)]];
%         yNose = [yNose,[yi(end-i);yi(end-i-1);yi(end-i-1);yi(end-i);yi(end-i)]];
%         zNose = [zNose,[j*0.0254/numMid;j*0.0254/numMid;(j+1)*0.0254/numMid;(j+1)*0.0254/numMid;j*0.0254/numMid]];
%         xNoseLattice = [xNoseLattice,[xi(end-i);xi(end-i-1);xi(end-i-1);xi(end-i)]];
%         yNoseLattice = [yNoseLattice,[yi(end-i);yi(end-i-1)-(yi(end-i-1)-yi(end-i))/4;yi(end-i-1)-(yi(end-i-1)-yi(end-i))/4;yi(end-i)]];
%         zNoseLattice = [zNoseLattice,[j*0.0254/numMid;j*0.0254/numMid;(j+1)*0.0254/numMid;(j+1)*0.0254/numMid]];
%         xNosem = [xNosem,xi(end-i)+c/4];
%         yNosem = [yNosem,(yi(end-i-1)-yi(end-i))/4+yi(end-i)];
%         zNosem = [zNosem,(j+.5)*0.0254/numMid];
%         sNose = [sNose,pi];
%     end
% end

[n,m] = size(xNose);
rotX = [1,0,0;0,cos(pi),-sin(pi);0,sin(pi),cos(pi)];
tempX = zeros(5,1);
tempY = zeros(5,1);
tempZ = zeros(5,1);
tempXL = zeros(4,1);
tempYL = zeros(4,1);
tempZL = zeros(4,1);

for i = 1:m
    for j = 1:n
        temp = rotX*[xNose(j,i);yNose(j,i);zNose(j,i)];
        tempX(j) = temp(1);
        tempY(j) = temp(2);
        tempZ(j) = temp(3);
    end
    xNose = [xNose,tempX];
    yNose = [yNose,tempY];
    zNose = [zNose,tempZ];
end

[n,m] = size(xNoseLattice);
for i = 1:m
    for j = 1:n
        tempL = rotX*[xNoseLattice(j,i);yNoseLattice(j,i);zNoseLattice(j,i)];
        tempXL(j) = tempL(1);
        tempYL(j) = tempL(2);
        tempZL(j) = tempL(3);
    end
    xNoseLattice = [xNoseLattice,tempXL];
    yNoseLattice = [yNoseLattice,tempYL];
    zNoseLattice = [zNoseLattice,tempZL];
end

tempM = rotX*[xNosem;yNosem;zNosem];
xNosem = [xNosem,tempM(1,:)];
yNosem = [yNosem,tempM(2,:)];
zNosem = [zNosem,tempM(3,:)];

sNose = [sNose,sNose];

zNose = zNose+zb*ones(size(zNose));
zNoseLattice = zNoseLattice+zb*ones(size(zNoseLattice));
zNosem = zNosem+zb*ones(size(zNosem));

yNose = yNose';
zNose = zNose';
xNosem = -(xNosem'-noseH);
yNosem = yNosem';
zNosem = zNosem';
xNoseLattice = -(xNoseLattice'-noseH);
yNoseLattice = yNoseLattice';
zNoseLattice = zNoseLattice';
xNose = -(xNose'-noseH);