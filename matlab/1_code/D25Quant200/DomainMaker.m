

D = 25; save('D.mat','D') %Diameter of the domain measured in units of Ro
quant = 200;save('quant.mat','quant')%minimal number of intervals covered by a radius

%geometry of Domain
nr = ceil(D*quant/2); save('nr.mat','nr')% number of intervals where integration in needed, which is also the number of points 
%along the radial axis where the function is possibly non-zero
dr = D/(2*nr); save('dr.mat','dr')%mesh resolution in units of Ro
nlmax = quant+1; save('nlmax.mat','nlmax')%max number of contact points possible for a non-deforming impactor
r = 0:dr:D/2; save('r.mat','r')%radial mesh
rn = 0:nr+1; save('rn.mat','rn')%radial mesh in units of dr
xplot = [-fliplr(r(2:nr+1)),r]; save('xplot.mat','xplot')%xaxis, used only for graphing

%Aux variables
zs = zeros(1,nlmax); %bottom of the ball/drop plus and infinite height away from the drop in units of drop radii
zs(1) = 0; 
for j=1:quant
    zs(j+1)=1-sqrt(1-j^2*dr^2);
end
zs=[zs,100*ones(1,nr-quant-1)]'; save('zs.mat','zs')

%Sistem matrixes
subdiag = [[zeros(1,nr-1);eye(nr-1)],zeros(nr,1)];
Dr = (-subdiag+subdiag')/2;
OneOverR = abs([0,r(2:nr).^(-1)]); save('OneOverR.mat','OneOverR')
DrOverR = zeros(nr,nr);
for i=2:nr
    DrOverR(i,:) = OneOverR(i)*Dr(i,:);
end

Derivr = Dr/dr; 
Derivr(1,:) = zeros(1,size(Derivr,2)); save('Derivr.mat','Derivr')
Derivrr = (subdiag+subdiag'-2*eye(nr))/dr^2; 
Derivrr(1,:) = [-2 2,zeros(1,nr-2)]/dr^2; save('Derivrr.mat','Derivrr')
Delta = Derivrr+DrOverR/dr;
Delta(1,:) = [-4 4,zeros(1,nr-2)]/dr^2;
save('Delta.mat','Delta')%Laplacian matrix

% tanDrop = zeros(1,nlmax); %tangent to the sphere at a given radial position
% tanDrop(1) = 0;
% for k = 2:nlmax
%     tanDrop(k) = dr*(k-1)/sqrt(1-(k-1)^2*dr^2);
% end
% angleDrop = atan(tanDrop); %angle of tangency of the drop at a given position
% save('angleDrop.mat','angleDrop')

tanDropMP = zeros(1,nlmax); %tangent to the sphere at a given radial position
for kk=1:nlmax-1
    tanDropMP(kk) = dr*(kk-.5)/sqrt(1-(kk-.5)^2*dr^2);
end
angleDropMP = atan(tanDropMP(1:nlmax-1)); %angle of tangency of the drop at a given position
angleDropMP = [angleDropMP,pi/2];
save('angleDropMP.mat','angleDropMP')

% Int = zeros(1,nlmax);
% Int(1) = 1/3;
% for k=2:nlmax
%     Int(k) = (2*(k-1));
% end
% Int = pi*dr^2*Int; %Integral functional, aproximated as lin by parts
% save('Int.mat','Int')

IntMat = zeros(2*nlmax,2*nlmax);
IntMat(1,1) = 1/12;
for ii=2:2*nlmax
    IntMat(ii,1) = 1/3; 
    for jj = 2:ii-1
        IntMat(ii,jj) = 2*jj-2;
    end
    IntMat(ii,ii) = 1.5*ii-21/12;
end
IntMat = pi*dr^2*IntMat;
save('IntMat.mat','IntMat')

%Drop geometry for graphing
xdrop = cos(2*pi*(0:.01:1)); %x coordinates of a sphere in dimensionless units, used to make movies
save('xdrop.mat','xdrop')
zdrop = sin(2*pi*(0:.01:1))+1; %z coordinates of the sphere in dimensionless units, used to make movies
save('zdrop.mat','zdrop')

% %Finding the DTN operator
% numer = ceil(pi*D/dr);
% if mod(numer,2)==1
%     numer = numer+1;
% end
% dtheta = 2*pi/numer;%use pi/even number
% %Integrating in the big circle
% DTN=zeros(nr,nr);
% for k=1:nr
%      k/nr
%     for i=1:rn(k)+nr
%         for l=dtheta/2:dtheta:pi-dtheta/4
%             radn = abs(sqrt((rn(k)+i*cos(l))^2+(i*sin(l))^2));
%             i1 = floor(radn);
%             w1 = min(max(0,radn - i1),1);
%             if i1 < .5
%                 DTN(k,i1+1) = DTN(k,i1+1) - (3*w1^3/4-7*w1^2/4+1)*    2*(1/(i-1/2)-1/(i+1/2));
%                 DTN(k,i1+2) = DTN(k,i1+2) - (-w1+2)*w1^2*             2*(1/(i-1/2)-1/(i+1/2));
%                 DTN(k,i1+3) = DTN(k,i1+3) - (w1-1)*w1^2/4*            2*(1/(i-1/2)-1/(i+1/2));
%             elseif i1<nr
%                 DTN(k,i1)   = DTN(k,i1)   - (-w1^2/6+w1/2-1/3)*w1*    2*(1/(i-1/2)-1/(i+1/2));
%                 DTN(k,i1+1) = DTN(k,i1+1) - (w1^3/2-w1^2-w1/2+1)*     2*(1/(i-1/2)-1/(i+1/2));
%                 if i1<nr-1
%                     DTN(k,i1+2) = DTN(k,i1+2) - (-w1^2/2+w1/2+1)*w1*  2*(1/(i-1/2)-1/(i+1/2));
%                     if i1<nr-2
%                         DTN(k,i1+3) = DTN(k,i1+3) - (w1^2-1)*w1/6*    2*(1/(i-1/2)-1/(i+1/2));
%                     end
%                 end
%             end
%         end
%     end
% end
% DTN = DTN*dtheta/(2*pi*dr);
% %integrating phi(x) exactly (the constant cancels here when integrating form zero to 2*\pi)
% for k=1:nr
%     DTN(k,k) = DTN(k,k) + 2/dr;
% end
% save(['DTN_Ro',num2str(10*Ro),'nlmax',num2str(nlmax),'D',num2str(D),'.mat'],'DTN')
   




