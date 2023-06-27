% This script will loop between all folders and compute the relevant
% metrics of the simulation.

% STEP 3: Actually run the simulations. 

files_folder = dir("**/etaOri.mat");
for ii = 1:length(files_folder)
    cd(files_folder(ii).folder);
    
    while ~isfile("dr.mat")
        cd ..
    end
    load("dr.mat");
    cd(files_folder(ii).folder);
    
    % Check if etaOri exists (the center of the bath)
    if isempty(dir("oscillation*.mat")) == false % || true
        
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
        contact_idx = find(north<=2,1);
        flight_idx = find(north((contact_idx+1):end)>=2,1);
        N = size(oscillation_amplitudes, 1);
        south = z - (1 + sum(oscillation_amplitudes .* ((-ones(N, 1)).^((1:N)')), 1));
        
        tImpact = (tvec(contact_idx)+tvec(contact_idx+1))/2; 
        Uo = (vz(contact_idx)+vz(contact_idx+1))/2;
        tend = (tvec(contact_idx+flight_idx)+tvec(contact_idx+flight_idx-1))/2;
        Uend = (vz(contact_idx+flight_idx)+vz(contact_idx+flight_idx-1))/2;
        tcont = tend-tImpact;
        CRref = -Uend/Uo;
        
        max_def = min(south); if max_def == -1; max_def = NaN; fprintf("Error on %s \n", pwd); end
        
        L = diff(etas)/dr;
        max_gradient = max(abs(L(:)));
        if max_gradient > 1; warning("Gradient too big for %s", pwd); end
        save('simulation_postprocessing.mat', "Uo", "tend", ...
            "Uend", "max_def", "CRref", "tcont", "max_gradient");

    end
end
