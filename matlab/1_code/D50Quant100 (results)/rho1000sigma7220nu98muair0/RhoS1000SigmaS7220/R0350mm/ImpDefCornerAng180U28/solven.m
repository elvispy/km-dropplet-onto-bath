function [etaprob,phiprob,zprob,vzprob,psprob,errortan] = solven(nl,dt,zo,vzo,etao,phio,nr,dr,Re,Delta,DTN,...
    Fr,We,Ma,zs,Int,angleDrop)
freeplaces = nl+1:nr;
pressedplaces = 1:nl;
b = [etao;phio];

Deltaprime = [zeros(nl,nr);Delta(nl+1:nr,:)];
Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
    [dt*(eye(nr)/Fr-Deltaprime/We),eye(nr)-dt*2*Delta/Re]];
bmod = b-Sist(:,pressedplaces)*zs(pressedplaces)...
        + [zeros(nr,1);2*dt*ones(nl,1)/We;zeros(nr-nl,1)];

ds = [[Sist(:,[freeplaces,nr+1:2*nr]),...
    [zeros(nr,nl);dt*eye(nl);zeros(nr-nl,nl)],...
    zeros(2*nr,1),Sist(:,pressedplaces)*ones(nl,1)];
    [zeros(1,2*nr-nl),-dt*Int(1:nl)/Ma,1,0];
    [zeros(1,2*nr-nl),-dt^2*Int(1:nl)/(2*Ma),0,1]]...
    \[bmod;vzo-dt/Fr;zo+vzo*dt-dt^2/(2*Fr)];

etaprob(freeplaces,1) = ds(1:size(freeplaces,2));
phiprob = ds(size(freeplaces,2)+1:2*nr-nl);
psprob = ds(2*nr-nl+1:2*nr);
vzprob = ds(2*nr+1);
zprob = ds(2*nr+2);
etaprob(pressedplaces,1) = zprob+zs(pressedplaces);

check = (etaprob(nl+1:end)>(zprob+zs(nl+1:end)));
if sum(check)>.5
    errortan = 4;
else
    taneffect = (etaprob(nl+1)-etaprob(nl))/dr;
    ataneffect = atan(taneffect);
    errortan = angleDrop(nl)-ataneffect;
end