% This script will loop between all folders and compute the relevant
% metrics of the simulation.

% STEP 3: Actually run the simulations. 

files_folder = dir("**/etaOri.mat");
for ii = 1:length(files_folder)
    cd(files_folder(ii).folder);
    
    % Check if etaOri exists (the center of the bath)
    if isempty(dir("oscillation*.mat")) == false
        
        try
            load('U0.mat');
        
        catch
            load("ProblemConditions.mat");
        end
        load('vz.mat'); 
        Vo = abs(vz(1));
        if exist('etas.m', 'file')
            load('etas.m');
        else
            files = dir(fullfile(pwd, "etaMatPer*.mat"));
            N = length(files);
            etaAux = [];
            for i = 1:N
                load(files(i).name);
                etaAux = [etaAux, etaMatPer];
            end
            etaMatPer = etaAux; etas = etaAux; save('etas.mat', 'etas');
        end
        load('z.mat')
        load('etaOri.mat')
        load('tvec.mat')
        load('oscillation_amplitudes.mat');
     
        
        north = z + 1 + sum(oscillation_amplitudes, 1);
        index1 = find(north<=2,1);
        index2 = find(north((index1+1):end)>=2,1);
        N = size(oscillation_amplitudes, 1);
        south = z - (1 + sum(oscillation_amplitudes .* ((-ones(N, 1)).^((1:N)')), 1));
        
        tImpact = (tvec(index1)+tvec(index1+1))/2; 
        Uo = (vz(index1)+vz(index1+1))/2;
        tend = (tvec(index1+index2-1)+tvec(index1+index2-2))/2;
        Uend = (vz(index1+index2-1)+vz(index1+index2-2))/2;
        tcont = tend-tImpact;
        CRref = -Uend/Uo;
        
        max_def = min(south); if max_def == -1; max_def = NaN; fprintf("Error on %s \n", pwd); end

        save('simulation_postprocessing.mat', "Uo", "tend", ...
            "Uend", "max_def", "CRref", "tcont");

    end
end
