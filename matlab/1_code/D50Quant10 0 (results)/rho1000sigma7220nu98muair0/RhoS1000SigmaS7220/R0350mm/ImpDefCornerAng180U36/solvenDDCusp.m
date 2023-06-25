function [etaprob,phiprob,zprob,vzprob,psprob,errortan] = solvenDDCusp(nlprev,nl,dt,zo,vzo,etao,phio,...
    nr,dr,Re,Delta,DTN,Fr,We,Ma,zs,Int,angleDrop,Cang,WeSB,Rv)

b = [etao;phio];
if nlprev < nl
    AngC = min(angleDrop(nl),pi-Cang);
else
    AngC = pi-Cang;
end

Deltaprime = [zeros(nl,nr);Delta(nl+1:nr,:)];
Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
    [dt*(eye(nr)/Fr-Deltaprime/We),eye(nr)-dt*2*Delta/Re]];
bmod = b-Sist(:,1:nl)*(zs(1:nl)+Rv);

Mat =  [[Sist(:,nl+1:2*nr),...
    [zeros(nr,nl);dt*eye(nl);zeros(nr-nl,nl)],...
    zeros(2*nr,1),Sist(:,1:nl)*ones(nl,1)];
    [zeros(1,2*nr-nl),-dt*Int(1:nl)/Ma,1,0];
    [zeros(1,2*nr-nl),-dt^2*Int(1:nl)/(2*Ma),0,1]];
% cond(Mat)
indep = [bmod;vzo-dt/Fr+dt*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(2*WeSB);...
    zo+vzo*dt-dt^2/(2*Fr)+dt^2*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(4*WeSB)];

ds = Mat\indep;
etaprob(nl+1:nr,1) = ds(1:nr-nl);
phiprob = ds(nr-nl+1:2*nr-nl);
psprob = ds(2*nr-nl+1:2*nr);
vzprob = ds(2*nr+1);
zprob = ds(2*nr+2);
etaprob(1:nl,1) = zprob+Rv+zs(1:nl);

check = (etaprob(nl+1:end)>(zprob+Rv+zs(nl+1:end)));
if sum(check)>.5
    errortan = 4;
else
    taneffect = (etaprob(nl+1)-etaprob(nl))/dr;
    ataneffect = atan(taneffect);
    errortan = angleDrop(nl)-ataneffect-AngC;
end