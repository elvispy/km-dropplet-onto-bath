function [etaprob,phiprob,zprob,vzprob,psprob,errortan] = getNextStep(nlprev,nl,dt,zo,vzo,etao,phio,...
    zs,Rv, PROBLEM_CONSTANTS, nlmax)

if nl < 0 || nl > nlmax
    etaprob = 0;
    phiprob = 0;
    zprob = 0;
    vzprob = 0;
    psprob = 0;
    errortan = Inf;
else
    We = PROBLEM_CONSTANTS.We; Delta = PROBLEM_CONSTANTS.Delta;
    Re = PROBLEM_CONSTANTS.Re; DTN = PROBLEM_CONSTANTS.DTN;
    nr = PROBLEM_CONSTANTS.nr; IntMat = PROBLEM_CONSTANTS.IntMat(max(nl, 1), :);
    angleDropMP = PROBLEM_CONSTANTS.angleDropMP; Cang = PROBLEM_CONSTANTS.Cang;
    Dr = PROBLEM_CONSTANTS.Dr; Ma = PROBLEM_CONSTANTS.Ma;
    dr = PROBLEM_CONSTANTS.dr; Fr = PROBLEM_CONSTANTS.Fr;

    if nl == 0
        % ??? this script seems wrong
        b = [etao;phio];
        
        Sist = [eye(nr)-dt*2*Delta/Re,-dt*DTN;...
               [dt*(eye(nr)/Fr-Delta/We),eye(nr)-dt*2*Delta/Re]];
        c = (Sist\b);
        
        etaprob = c(1:nr);
        phiprob = c(nr+1:2*nr);
        vzprob = vzo-dt/Fr;
        
        %Finding new eta
        zprob = zo+vzo*dt-dt^2/(2*Fr);
        check = (etaprob>(zprob+zs+Rv));
        errortan = 0;
        if sum(check)>.5
            errortan = Inf;
        end
        psprob = [];
    else % There are contact points
        [etaprob,phiprob,zprob,vzprob,psprob,errortan] = solvenDDCusp(nlprev,nl,dt,zo,vzo,etao,phio,...
            nr,dr,Re,Delta,DTN,Fr,We,Ma,zs,IntMat,angleDropMP,Cang,Dr,Rv);
    end
end
