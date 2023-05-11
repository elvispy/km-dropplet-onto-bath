

tic
maxtime = 10*60*60;
poolnum = 4; %Size of the parallel pool
workerload = 10;%Number of tasks that we expect worker to run before saving
%load('runNumber.mat','runNumber')
runNumber = 0;
if runNumber == 0
    
    load('nr.mat','nr')
    load('D.mat','D')
    load('rn.mat','rn')
    load('dr.mat','dr')
    
    refp = 10;
    save('refp.mat','refp')

    %Finding the DTN operator
    drp = dr/refp;
    numer = ceil(pi*D/drp);
    if mod(numer,2)==1
        numer = numer+1;
    end
    dtheta = 2*pi/numer;%use pi/even number
    DTNnew345=zeros(nr,nr);

    
    %% Integrating away from the singularity
    k = 1;
    for i=2:(rn(k)+nr+1)
        i1 = round(i);
        if i1<nr+1
            DTNnew345(k,i1)   = DTNnew345(k,i1)           - (( - i^2/2 -  i -1/3)*log((i+1)/(i)) +( i^2/6 +i/2 +1/3)     *(1-i/(i+1)) -(i+.5)/6 +i/2 +1/2);
        end
        if i1<nr
            DTNnew345(k,i1+1) = DTNnew345(k,i1+1)         - (( 3*i^2/2 +2*i -1/2)*log((i+1)/(i)) +(-i^2/2 -i   +1/2 +1/i)*(1-i/(i+1)) +(i+.5)/2 -3*i/2   -1);
            if i1<nr-1
                DTNnew345(k,i1+2) = DTNnew345(k,i1+2)     - ((-3*i^2/2 -  i +1  )*log((i+1)/(i)) +( i^2/2 +i/2 - 1)      *(1-i/(i+1)) -(i+.5)/2 +3*i/2   +1/2);
                if i1<nr-2
                    DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - ((   i^2/2      -1/6)*log((i+1)/(i)) +(-i^2        +1)/6     *(1-i/(i+1)) +(i+.5)/6 -i/2);
                end
            end
        end
    end
    DTNnew345(1,1) = DTNnew345(1,1) + 1/2;
    DTNnew345(1,:) = DTNnew345(1,:)/dr;
    %Integrating the vincinity of the origin (= sing)
    DTNnew345(1,1)   = DTNnew345(k,1) +209/(54*dr);%- (-15+    4)/(6*dr);
    DTNnew345(1,2)   = DTNnew345(k,2) -29/(6*dr);%+ (-16+4/3*4)/(6*dr);
    DTNnew345(1,3)   = DTNnew345(k,3) + 7/(6*dr);%+ (  1-1/3*4)/(6*dr);
    DTNnew345(1,4)   = DTNnew345(k,4) - 11/(54*dr);%- (-15+    4)/(6*dr);

    
    %% Integrating the vincinity of the origin
    k = 2;
    for i = 1:2*refp
        Kern = 2*(1/(i-1/2)-1/(i+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
            x1 = i*cos(l)/refp;
            posr = radn-rn(k);
            DTNnew345(k,k)   = DTNnew345(k,k) - ((  -6- 3*posr+9*posr^2+ 3*posr^3- 3*posr^4)*posr+6*x1)/72   *Kern;
            DTNnew345(k,k-1) = DTNnew345(k,k-1) - ((-88+48*posr+62*posr^2-12*posr^3-10*posr^4)*posr+88*x1)/72   *Kern;
            DTNnew345(k,k)   = DTNnew345(k,k)   - ((  72-90*posr-90*posr^2+18*posr^3+18*posr^4)*posr- 72*x1)/72    *Kern;
            if k<nr-.5
                DTNnew345(k,k+1) = DTNnew345(k,k+1) - (( 24+48*posr+18*posr^2-12*posr^3-6*posr^4)*posr-24*x1)/72   *Kern;
                if k < nr-1.5
                    DTNnew345(k,k+2) = DTNnew345(k,k+2) - ((  -2- 3*posr+ posr^2+ 3*posr^3+ posr^4)*posr+ 2*x1)/72   *Kern;
                end
            end    
        end
    end
    %% Integrating away from the singularity
    for i=2*refp+1:(rn(k)+nr)*refp
        Kern = 2*(1/(i-1/2)-1/(i+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
            i1 = floor(radn);
            w1 = min(max(0,radn - i1),1);
            if i1 < .5
                    DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                    DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1+2)*w1^2*             Kern;
                    DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1-1)*w1^2/4*            Kern;
            elseif i1<nr
                DTNnew345(k,i1)   = DTNnew345(k,i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                if i1<nr-1
                    DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                    if i1<nr-2
                        DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1^2-1)*w1/6*    Kern;
                    end
                end
            end
        end
    end
    DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
    DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp);   


    %% Integrating the vincinity of the origin
    k = 3;
    for i = 1:2*refp
        Kern = 2*(1/(i-1/2)-1/(i+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
            x1 = i*cos(l)/refp;
            posr = radn-rn(k);
            DTNnew345(k,k-2) = DTNnew345(k,k-2) - (( 124- 12*posr-149*posr^2+ 12*posr^3+ 25*posr^4)*posr-124*x1)/288   *Kern;
            DTNnew345(k,k-1) = DTNnew345(k,k-1) - ((-384+192*posr+288*posr^2- 48*posr^3- 48*posr^4)*posr+384*x1)/288   *Kern;
            DTNnew345(k,k)   = DTNnew345(k,k)   + ((-  4+ 10*posr+  5*posr^2-  2*posr^3-    posr^4)*posr+  4*x1)/8     *Kern;
            if k<nr-.5
                DTNnew345(k,k+1) = DTNnew345(k,k+1) - (( 128+192*posr+ 32*posr^2- 48*posr^3- 16*posr^4)*posr-128*x1)/288   *Kern;
                if k < nr-1.5
                    DTNnew345(k,k+2) = DTNnew345(k,k+2) - ((- 12- 12*posr+  9*posr^2+ 12*posr^3+ 3*posr^4)*posr+ 12*x1)/288   *Kern;
                end
            end    
        end
    end
    %% Integrating away from the singularity
    for i=2*refp+1:(rn(k)+nr)*refp
        Kern = 2*(1/(i-1/2)-1/(i+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
            i1 = floor(radn);
            w1 = min(max(0,radn - i1),1);
            if i1 < .5
                    DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                    DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1+2)*w1^2*             Kern;
                    DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1-1)*w1^2/4*            Kern;
            elseif i1<nr
                DTNnew345(k,i1)   = DTNnew345(k,i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                if i1<nr-1
                    DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                    if i1<nr-2
                        DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1^2-1)*w1/6*    Kern;
                    end
                end
            end
        end
    end
    DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
    DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp); 
    
    % Shut down existing parallel pool.
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    else
        disp('No pool exists.')
    end
    % Create a new parallel pool.
    if isempty(gcp('nocreate'))
        parpool(poolnum)
    else
        disp('A pool already exists.')
    end
    for ll = 4:poolnum*workerload:nr 
        parfor k = ll:min(ll+poolnum*workerload-1,nr)
            k/nr
            Line = zeros(1,nr);
            for i = 1:2*refp
                Kern = 2*(1/(i-1/2)-1/(i+1/2));
                for l=dtheta/2:dtheta:pi-dtheta/4
                    radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
                    x1 = i*cos(l)/refp;
                    posr = radn-rn(k);
                    Line(k-2) = Line(k-2) - (  ( 2-  posr-2*posr^2+  posr^3)*posr- 2*x1)/24   *Kern;
                    Line(k-1) = Line(k-1) - (4*(-4+4*posr+  posr^2-  posr^3)*posr+16*x1)/24   *Kern;
                    Line(k)   = Line(k)   - (posr^2-5)*posr^2/4                               *Kern;
                    if k<nr-.5
                        Line(k+1) = Line(k+1) -  (4*(4+4*posr- posr^2- posr^3)*posr-16*x1)/24 *Kern;
                        if k < nr-1.5
                            Line(k+2) = Line(k+2) - ((-2-posr+2*posr^2+posr^3)*posr+ 2*x1)/24 *Kern;
                        end
                    end    
                end
            end     
            for i=2*refp+1:(rn(k)+nr)*refp
                Kern = 2*(1/(i-1/2)-1/(i+1/2));
                for l=dtheta/2:dtheta:pi-dtheta/4
                    radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
                    i1 = floor(radn);
                    w1 = min(max(0,radn - i1),1);
                    if i1 < .5
                        Line(i1+1) = Line(i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                        Line(i1+2) = Line(i1+2) - (-w1+2)*w1^2*             Kern;
                        Line(i1+3) = Line(i1+3) - (w1-1)*w1^2/4*            Kern;
                    elseif i1<nr
                        Line(i1)   = Line(i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                        Line(i1+1) = Line(i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                        if i1<nr-1
                            Line(i1+2) = Line(i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                            if i1<nr-2
                                Line(i1+3) = Line(i1+3) - (w1^2-1)*w1/6*    Kern;
                            end
                        end
                    end
                end
            end
            Line(:) = dtheta/(2*pi*drp)*Line(:);
            Line(k) = Line(k) + 2/(4*dr+drp);

            DTNnew345(k,:) = Line;
        end
        
        
        save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        
        runtime = toc;
        if runtime >= maxtime
            reStartAt = ll + poolnum*workerload;
            save('reStartAt.mat','reStartAt')
            runNumber = 1;
            save('runNumber.mat','runNumber')
            % Shut down the pool.
            if ~isempty(gcp('nocreate'))
                delete(gcp);
            else
                disp('No pool exists.')
            end
            % exit
        end 
    end
    save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
    reStartAt = nr+1;
    save('reStartAt.mat','reStartAt')
    runNumber = 1;
    save('runNumber.mat','runNumber')
    % Shut down the pool.
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    else
        disp('No pool exists.')
    end  
    % exit
elseif runNumber > 0
    load('nr.mat','nr')
    load('D.mat','D')
    load('rn.mat','rn')
    load('dr.mat','dr')
    load('refp.mat','refp')
    load('reStartAt.mat','reStartAt')
    
    drp = dr/refp;
    numer = ceil(pi*D/drp);
    if mod(numer,2)==1;
        numer = numer+1;
    end
    dtheta = 2*pi/numer;%use pi/even number
    
    if reStartAt > nr
        % exit
    else
        load(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        for ll = reStartAt:poolnum*workerload:nr 
            parfor k = ll:min(ll+poolnum*workerload-1,nr)
                k/nr
                Line = zeros(1,nr);
                for i = 1:2*refp
                    Kern = 2*(1/(i-1/2)-1/(i+1/2));
                    for l=dtheta/2:dtheta:pi-dtheta/4
                        radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
                        x1 = i*cos(l)/refp;
                        posr = radn-rn(k);
                        Line(k-2) = Line(k-2) - (  ( 2-  posr-2*posr^2+  posr^3)*posr- 2*x1)/24   *Kern;
                        Line(k-1) = Line(k-1) - (4*(-4+4*posr+  posr^2-  posr^3)*posr+16*x1)/24   *Kern;
                        Line(k)   = Line(k)   - (posr^2-5)*posr^2/4                               *Kern;
                        if k<nr-.5
                            Line(k+1) = Line(k+1) -  (4*(4+4*posr- posr^2- posr^3)*posr-16*x1)/24 *Kern;
                            if k < nr-1.5
                                Line(k+2) = Line(k+2) - ((-2-posr+2*posr^2+posr^3)*posr+ 2*x1)/24 *Kern;
                            end
                        end    
                    end
                end     
                for i=2*refp+1:(rn(k)+nr)*refp
                    Kern = 2*(1/(i-1/2)-1/(i+1/2));
                    for l=dtheta/2:dtheta:pi-dtheta/4
                        radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
                        i1 = floor(radn);
                        w1 = min(max(0,radn - i1),1);
                        if i1 < .5
                            Line(i1+1) = Line(i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                            Line(i1+2) = Line(i1+2) - (-w1+2)*w1^2*             Kern;
                            Line(i1+3) = Line(i1+3) - (w1-1)*w1^2/4*            Kern;
                        elseif i1<nr
                            Line(i1)   = Line(i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                            Line(i1+1) = Line(i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                            if i1<nr-1
                                Line(i1+2) = Line(i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                                if i1<nr-2
                                    Line(i1+3) = Line(i1+3) - (w1^2-1)*w1/6*    Kern;
                                end
                            end
                        end
                    end
                end
                Line(:) = dtheta/(2*pi*drp)*Line(:);
                Line(k) = Line(k) + 2/(4*dr+drp);
                DTNnew345(k,:) = Line;
            end
            save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
            runtime = toc;
            if runtime >= maxtime
                reStartAt = ll + poolnum*workerload;
                save('reStartAt.mat','reStartAt')
                runNumber = runNumber+1;
                save('runNumber.mat','runNumber')
                % Shut down the pool.
                if ~isempty(gcp('nocreate'))
                    delete(gcp);
                else
                    disp('No pool exists.')
                end
                % exit
            end
        end
        save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        reStartAt = nr+1;
        save('reStartAt.mat','reStartAt')
        runNumber = runNumber+1;
        save('runNumber.mat','runNumber')
        % Shut down the pool.
        if ~isempty(gcp('nocreate'))
            delete(gcp);
        else
            disp('No pool exists.')
        end  
        % exit
    end
end
