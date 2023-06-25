function [etaprob0,phiprob0,zprob0,vzprob0,errortan]=solve0(dt,z,vz,etao,phio,nr,Re,Delta,DTN,Fr,We,zs)

b = [etao;phio];

Sist = [eye(nr)-dt*2*Delta/Re,-dt*DTN;...
    [dt*(eye(nr)/Fr-Delta/We),eye(nr)-dt*2*Delta/Re]];
c = (Sist\b);

etaprob0 = c(1:nr);
phiprob0 = c(nr+1:2*nr);

vzprob0 = vz-dt/Fr;
%Finding new eta
zprob0 = z+vz*dt-dt^2/(2*Fr);
check = (etaprob0>(zprob0+zs));
errortan = 0;
if sum(check)>.5
    errortan = 4;
end