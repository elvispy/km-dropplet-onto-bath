function [etaprob,phiprob,zprob,vzprob,psprob,errortan] = solvenDDCusp1(nlprev,nl,dt,zo,vzo,etao,phio,...
    nr,dr,Re,Delta,DTN,Fr,We,Ma,zs,Int,angleDrop,Cang,WeS,Rv)
freeplaces = nl+1:nr;
pressedplaces = 1:nl;
b = [etao;phio];
% if nlprev < nl
%     AngC = min(angleDrop(nl),pi-Cang);
% else
    AngC = pi-Cang;
% end

Deltaprime = [zeros(nl,nr);Delta(nl+1:nr,:)];
Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
    [dt*(eye(nr)/Fr-Deltaprime/We),eye(nr)-dt*2*Delta/Re]];
bmod = b-Sist(:,pressedplaces)*(zs(pressedplaces)+Rv);

Mat =  [[Sist(:,[freeplaces,nr+1:2*nr]),...
    [zeros(nr,nl);dt*eye(nl);zeros(nr-nl,nl)],...
    zeros(2*nr,1),Sist(:,pressedplaces)*ones(nl,1)];
    [zeros(1,2*nr-nl),-dt*Int(1:nl)/Ma,1,0];
    [zeros(1,2*nr-nl),-dt^2*Int(1:nl)/(2*Ma),0,1]];%     [zeros(1,2*nr),-dt,1]]...
Mat = [Mat(nl+1:nr,:);Mat(1:nl,:);Mat(nr+nl+1:2*nr,:);Mat(nr+1:nr+nl,:);Mat(end-1:end,:)];
rcond(Mat)
indep = [bmod;vzo-dt/Fr+dt*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(2*WeS);...%I need to change WeS when the two fluids have different surface tensions
    zo+vzo*dt-dt^2/(2*Fr)+dt^2*3*(nl-.5)*dr*sin(angleDrop(nl)-AngC)/(4*WeS)];
%     zo];
indep = [indep(nl+1:nr);indep(1:nl);indep(nr+nl+1:2*nr);indep(nr+1:nr+nl);indep(end-1:end)];
ds = Mat\indep;
etaprob(freeplaces,1) = ds(1:size(freeplaces,2));
phiprob = ds(size(freeplaces,2)+1:2*nr-nl);
psprob = ds(2*nr-nl+1:2*nr);
vzprob = ds(2*nr+1);
zprob = ds(2*nr+2);
etaprob(pressedplaces,1) = zprob+Rv+zs(pressedplaces);

check = (etaprob(nl+1:end)>(zprob+Rv+zs(nl+1:end)));
if sum(check)>.5
    errortan = 4;
else
    taneffect = (etaprob(nl+1)-etaprob(nl))/dr;
    ataneffect = atan(taneffect);
    errortan = angleDrop(nl)-ataneffect-AngC;
end