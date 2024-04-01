files = dir("**/*.mat");

for ii = 1:length(files)
    cd(files(ii).folder)
    a = load(files(ii).name);
    
    if isfield(a, 'var')
        extract_varname = regexp(files(ii).name,'(errored_)?(.*)\.mat','tokens','once');
        extract_varname = extract_varname{end};
        a.(extract_varname) = a.var;
        a = rmfield(a, 'var');
        save(files(ii).name, '-struct', 'a');
    end
    
end