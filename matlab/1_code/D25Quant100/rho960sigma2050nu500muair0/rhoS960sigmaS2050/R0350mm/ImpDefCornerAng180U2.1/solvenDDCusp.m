function [etaprob,phiprob,zprob,vzprob,psprob,errortan] = solvenDDCusp(nlprev,nl,dt,zo,vzo,etao,phio,...
    nr,dr,Re,Delta,DTN,Fr,We,~,zs,Int,angleDrop,Cang,Dr,Rv)
Ma = 4*pi/3; We = We/Dr;
b = [etao;phio];
if nlprev < nl
    AngC = min(angleDrop(nl),pi-Cang);
else
    AngC = pi-Cang;
end

%Deltaprime = [zeros(nl,nr);Delta(nl+1:nr,:)];

% Preparing the matrix (2x2 block)
Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
    [dt*(eye(nr)/Fr-Delta/We),eye(nr)-dt*2*Delta/Re]];
bmod = b-Sist(:,1:nl)*(zs(1:nl)+Rv);

% Completing the system
Mat =  [[Sist(:,nl+1:2*nr),...
    [zeros(nr,nl);dt*eye(nl)*Dr;zeros(nr-nl,nl)],...
    zeros(2*nr,1),Sist(:,1:nl)*ones(nl,1)];
    [zeros(1,2*nr-nl),-dt*Int(1:nl)/Ma,1 ,0];
    [zeros(1,2*nr-nl),-zeros(1, nl)  ,-dt,1]];
% cond(Mat)
%indep = [bmod;vzo-dt/Fr+dt*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(2*We);...
%    zo+vzo*dt-dt^2/(2*Fr)+dt^2*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(4*We)];

indep = [bmod;vzo-dt/Fr;zo];

% Solving 
ds = Mat\indep;
etaprob(nl+1:nr,1) = ds(1:nr-nl);
phiprob = ds(nr-nl+1:2*nr-nl);
psprob = ds(2*nr-nl+1:2*nr);
vzprob = ds(end-1);
zprob = ds(end);
etaprob(1:nl,1) = zprob+Rv+zs(1:nl);

check = (etaprob(nl+1:end)>(zprob+Rv+zs(nl+1:end)));
if sum(check)>.5
    errortan = 4;
else
    taneffect = (etaprob(nl+1)-etaprob(nl))/dr;
    ataneffect = atan(taneffect);
    errortan = angleDrop(nl)-ataneffect-AngC;
end