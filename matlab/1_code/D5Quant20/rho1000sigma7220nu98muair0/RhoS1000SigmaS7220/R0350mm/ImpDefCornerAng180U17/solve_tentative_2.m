function [etaprob, phiprob, zprob,vzprob,errortan, ...
            numlTent, etaTent, phiTent, psNew, zTent, vzTent, ...
            tvec, tentative_index, reduc] ...
            = solve_tentative_2(numl, tentative_index, dt, z, vz, etao, phio, nr, Re, Delta, DTN, Fr, We, zs, RvTent, ...
            errortan, etaprob, phiprob, zprob, vzprob, psprob, dr, Ma, IntMat, angleDropMP, Cang, Dr, tvec, reduc, nlmaxTent, ...
            numlTent, etaTent, phiTent, psNew, zTent, vzTent)
        if numl(tentative_index) < .5 %i.e. if previously in flight (I need to define this as integer)
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),errortan(3,tentative_index+1)] = ...
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(3,tentative_index+1))<.5
                numlTent = 0;
                etaTent = etaprob(:,3);
                phiTent = phiprob(:,3);
                psNew = zeros(nlmaxTent,1);
                zTent = zprob(3);
                vzTent = vzprob(3);
            else
                co = find(numl(tentative_index:-1:1)~=1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1,4),errortan(4,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=2,1);
                [~,~,~,~,~,errortan(5,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                    numlTent = 1;
                    etaTent = etaprob(:,4);
                    phiTent = phiprob(:,4);
                    psNew = psprob(:,4);
                    zTent = zprob(4);
                    vzTent = vzprob(4);
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1;
                    reduc = 1;
                end
            end
        elseif numl(tentative_index)>.5 && numl(tentative_index)<1.5 % i.e. the last number of contact points was 1
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),errortan(2,tentative_index+1)] = ...
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(2,tentative_index+1))<.5
                numlTent = 0;
                etaTent = etaprob(:,2);
                phiTent = phiprob(:,2);
                psNew = zeros(nlmaxTent,1);
                zTent = zprob(2);
                vzTent = vzprob(2);
            else
                co = find(numl(tentative_index:-1:1)~=1,1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1,3),errortan(3,tentative_index+1)] = ...                     
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=2,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:2,4),errortan(4,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                    numlTent = 1;
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(tentative_index:-1:1)~=3,1);
                    [~,~,~,~,~,errortan(5,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                        numlTent = 2;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                end
            end
        elseif numl(tentative_index) > 1.5 && numl(tentative_index) < 2.5 %i.e. the last contact had two points
            [~,~,~,~,errortan(1,tentative_index+1)] = ...
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(1,tentative_index+1))<.5
                tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                tentative_index = tentative_index-1; 
                reduc = 1;
            else
                co = find(numl(tentative_index:-1:1)~=2,1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:2,3),errortan(3,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                    numlTent = 1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    co = find(numl(tentative_index:-1:1)~=3,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:3,4),errortan(4,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                        numlTent = 2;
                        etaTent = etaprob(:,3);
                        phiTent = phiprob(:,3);
                        psNew = psprob(:,3);
                        zTent = zprob(3);
                        vzTent = vzprob(3);
                    else
                        co = find(numl(tentative_index:-1:1)~=4,1);
                        [~,~,~,~,~,errortan(5,tentative_index+1)] = ...    
                            solvenDDCusp(numl(tentative_index-co+1),4,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                            We,Ma,zs,IntMat(4,:),angleDropMP,Cang,Dr,RvTent);
                        if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                            numlTent = 3;
                            etaTent = etaprob(:,4);
                            phiTent = phiprob(:,4);
                            psNew = psprob(:,4);
                            zTent = zprob(4);
                            vzTent = vzprob(4);
                        else
                            tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                            tentative_index = tentative_index-1; 
                            reduc = 1;
                        end
                    end
                end
            end
        elseif numl(tentative_index)>2.5 && numl(tentative_index)<nlmaxTent-1.5 %if the last number of contact points was far from the boundaries
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                    numlTent = numl(tentative_index)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                end
            else
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(tentative_index)+1,4),errortan(4,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)+1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(3,tentative_index+1))<abs(errortan(4,tentative_index+1))
                    numlTent = numl(tentative_index);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+2,1);%I think I don't need this and I can just replace the first argument of solven by numl(jj)
                    [~,~,~,~,~,errortan(5,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)+2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                        numlTent = numl(tentative_index)+1;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                end
            end
        elseif numl(tentative_index) > nlmax(tentative_index)-1.5 && numl(tentative_index) < nlmaxTent-.5 %i.e. if last number of contacted points is nlmax-1
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1))<abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                    numlTent = numl(tentative_index)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                end
            else
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(tentative_index)+1,4),errortan(4,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)+1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                    numlTent = numl(tentative_index);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    numlTent = numl(tentative_index)+1;
                    etaTent = etaprob(:,4);
                    phiTent = phiprob(:,4);
                    psNew = psprob(:,4);
                    zTent = zprob(4);
                    vzTent = vzprob(4);
                end
            end
        elseif numl(tentative_index) == nlmaxTent %i.e. if last number of contact points was nlmax
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                    numlTent = numl(tentative_index)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                end
            else
                numlTent = numl(tentative_index);
                etaTent = etaprob(:,3);
                phiTent = phiprob(:,3);
                psNew = psprob(:,3);
                zTent = zprob(3);
                vzTent = vzprob(3);
            end
        else
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < 4
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                    numlTent = numl(tentative_index)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                end
            else
                tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                tentative_index = tentative_index-1; 
                reduc = 1;
            end
        end        
        
end