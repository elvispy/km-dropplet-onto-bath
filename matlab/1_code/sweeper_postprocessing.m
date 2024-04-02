% This script will loop between all folders and compute the relevant
% metrics of the simulation.

% STEP 3: Actually run the simulations. 
diary ../0_data/manual/Logger/sweeper_postprocessing_logger.txt
disp("-------");
fprintf("%s \n %s", datestr(datetime()), mfilename('fullpath'));
files_folder = dir("**/etaOri.mat");
for ii = 1:length(files_folder)
    cd(files_folder(ii).folder);
    
    while ~isfile("dr.mat")
        cd ..
    end
    load("dr.mat");
    cd(files_folder(ii).folder);
    mypwd = split(pwd, "1_code"); mypwd = mypwd{2};
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
            load('etas.mat', 'etas');
            files = dir(fullfile(pwd, "etaMatPer*.mat"));
            N = length(files);
            if N > 0
                etaAux = [];
                for i = 1:N
                    load(files(i).name);
                    etaAux = [etaAux, etaMatPer];
                end
                etaMatPer = etaAux;
            else
                etaMatPer = etas;
            end
        end
        load("numl.mat");
        load('z.mat')
        load('etaOri.mat')
        load('tvec.mat')
        load('oscillation_amplitudes.mat');
     
        % COM-Based parameters
        north = z + 1; % Changed: Now is COM-based. OLD used: + sum(oscillation_amplitudes, 1); 
        
        contact_idx = find(north<=2,1);

        flight_idx = find(north((contact_idx+1):end)>2,1);
        N = size(oscillation_amplitudes, 1);
        south = z - (1 + sum(oscillation_amplitudes .* ((-ones(N, 1)).^((1:N)')), 1));
        
        % Adimensional parameters
        tImpact = (tvec(contact_idx)+tvec(contact_idx+1))/2; 
        Uo = (vz(contact_idx)+vz(contact_idx+1))/2;
        tend = (tvec(contact_idx+flight_idx)+tvec(contact_idx+flight_idx-1))/2;
        Uend = (vz(contact_idx+flight_idx)+vz(contact_idx+flight_idx-1))/2;
        tcont = tend-tImpact;
        CRref = -Uend/Uo;
        
        max_def = min(south); if max_def == -1; max_def = NaN; fprintf("Error on %s \n", mypwd); end
        if norm(oscillation_amplitudes) < 1e-5
           save("delete_me.m", "oscillation_amplitudes");
           disp("lol");
        end
        
        % Number-of-contact-points based parameters
        contact_indicator = (numl ~= 0);
        transition_indicator = diff(contact_indicator);
        contacts = (transition_indicator == 1);
        liftoffs = (transition_indicator == -1);
        %nnew = diff(numl);
        %liftoff = [nnew, NaN] == numl;
        %contact = [NaN, nnew] == -numl;
        idx_impact_theory = find(contacts);
        idx_end_theory    = find(liftoffs);
        if isempty(idx_end_theory)
            fprintf("No liftoff in %s\n", mypwd);
            liftoff_times = nan;
            t_cont_theory = nan;
            contact_times = nan;
            Uend_theory = nan;
            Uo_theory = nan;
        else
        
            Uo_theory = vz(idx_impact_theory(1));
            Uend_theory = vz(idx_end_theory(1));
            t_cont_theory = tvec(idx_end_theory(1)) - tvec(idx_impact_theory(1));
            contact_times = tvec(idx_impact_theory);
            liftoff_times = tvec(idx_end_theory);
        end
        
        if length(idx_end_theory) > 1
            fprintf("Second bounce in %s\n", mypwd);     
        end
        try    
            CR_theory = -Uend_theory/Uo_theory;
        catch
            CR_theory = NaN;
        end
        
        L = diff(etas)/dr;
        max_gradient = max(abs(L(:)));
        if max_gradient > 1; warning("Big gradient in %s\n", mypwd); end
        save('simulation_postprocessing.mat', "Uo", "tend", ...
            "Uend", "max_def", "CRref", "tcont", "max_gradient", "Uo_theory", ...
            'Uend_theory', 't_cont_theory', 'contact_times', 'liftoff_times', ...
            'CR_theory');

    end
end
diary off
