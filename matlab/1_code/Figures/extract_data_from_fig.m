function summary = extract_data_from_fig(figFile, outDir)
%EXTRACT_DATA_FROM_FIG  Extract X/Y(/Z) from real Line & Scatter series in a .fig.
%   SUMMARY = EXTRACT_DATA_FROM_FIG(FIGFILE)
%   SUMMARY = EXTRACT_DATA_FROM_FIG(FIGFILE, OUTDIR)
%
%   Improvements vs. naive traversal:
%     - uses findobj (no hidden handles)
%     - ignores legend/colorbar/utility axes
%     - keeps only legend-eligible series (IconDisplayStyle='on')
%     - drops "series" with < 2 finite points (kills legend proxies/caps)
%     - skips errorbar hggroup child lines
%     - de-duplicates identical (X,Y[,Z]) series per axes
%     - writes one CSV per real series, plus a MAT with everything
%
%   Returns a SUMMARY struct with counts and file paths.

    narginchk(1,2);
    if ~isfile(figFile)
        error('File not found: %s', figFile);
    end

    [fp, base, ~] = fileparts(figFile);
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(fp, [base '_extracted']);
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Open figure invisibly
    fig = openfig(figFile, 'invisible');
    cleaner = onCleanup(@() safeClose(fig));

    % Axes: visible data axes only
    axs = findobj(fig, 'Type','axes', ...
        '-not','Tag','legend', '-not','Tag','Legend', ...
        '-not','Tag','Colorbar', '-not','Tag','scribeOverlay');
    axs = axs(:).'; % row

    sanitize = @(s) regexprep(strtrim(string(s)), '[^A-Za-z0-9_\-]+', '_');
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    extracted = struct();
    extracted.sourceFig   = which(figFile);
    extracted.createdAt   = datetime('now');
    extracted.axes        = [];
    savedFiles = {};

    axCounter = 0;
    for ax = axs
        axCounter = axCounter + 1;
        axInfo = struct();
        axInfo.index = axCounter;
        axInfo.title = tryAxTitle(ax);

        % =======================
        % Filtered LINE objects
        % =======================
        rawLines = findobj(ax, 'Type','line');   % no hidden handles
        keepLines = gobjects(0);
        for h = reshape(rawLines,1,[])
            if ~isLegendOn(h), continue; end
            if isErrorbarChild(h), continue; end

            X = get(h,'XData'); Y = get(h,'YData');
            if ~hasAtLeastTwoFinitePoints(X, Y), continue; end

            keepLines(end+1) = h; %#ok<AGROW>
        end
        keepLines = dedupByData(keepLines, 'line');

        % Export lines
        axInfo.lines = [];
        for j = 1:numel(keepLines)
            h = keepLines(j);
            [X, Y, Z] = getXYZ(h);
            dname = sanitize(tryGet(h,'DisplayName',sprintf('line_%d', j)));
            T = makeLineTable(X, Y, Z);
            if isempty(T), continue; end

            baseName = sprintf('%s_axes%02d_line%02d_%s_%s', ...
                base, axCounter, j, char(dname), timestamp);
            csvPath = fullfile(outDir, [baseName '.csv']);
            writetable(T, csvPath);
            savedFiles{end+1} = csvPath; %#ok<AGROW>

            s = struct('index', j, 'displayName', string(dname), ...
                       'X', X(:), 'Y', Y(:), 'Z', Z(:), 'fileCSV', csvPath);
            axInfo.lines = [axInfo.lines; s]; %#ok<AGROW>
        end

        % ==========================
        % Filtered SCATTER objects
        % ==========================
        rawSc = findobj(ax, 'Type','scatter');   % includes scatter & scatter3
        keepSc = gobjects(0);
        for h = reshape(rawSc,1,[])
            if ~isLegendOn(h), continue; end
            X = get(h,'XData'); Y = get(h,'YData');
            if ~hasAtLeastTwoFinitePoints(X, Y), continue; end
            keepSc(end+1) = h; %#ok<AGROW>
        end
        keepSc = dedupByData(keepSc, 'scatter');

        % Export scatter
        axInfo.scatter = [];
        for j = 1:numel(keepSc)
            h = keepSc(j);
            [X, Y, Z] = getXYZ(h);
            S = tryGet(h,'SizeData',[]);
            C = tryGet(h,'CData',[]);
            [C1,C2,C3] = normalizeC(C, numel(X));
            S = normalizeSize(S, numel(X));

            dname = sanitize(tryGet(h,'DisplayName',sprintf('scatter_%d', j)));

            T = table();
            if ~isempty(X), T.X = X(:); end
            if ~isempty(Y), T.Y = Y(:); end
            if ~isempty(Z), T.Z = Z(:); end
            if ~isempty(S), T.Size = S(:); end
            if ~isempty(C1), T.C1 = C1(:); end
            if ~isempty(C2), T.C2 = C2(:); end
            if ~isempty(C3), T.C3 = C3(:); end
            if isempty(T), continue; end

            baseName = sprintf('%s_axes%02d_scatter%02d_%s_%s', ...
                base, axCounter, j, char(dname), timestamp);
            csvPath = fullfile(outDir, [baseName '.csv']);
            writetable(T, csvPath);
            savedFiles{end+1} = csvPath; %#ok<AGROW>

            s = struct('index', j, 'displayName', string(dname), ...
                       'X', X(:), 'Y', Y(:), 'Z', Z(:), ...
                       'SizeData', S, 'CData', struct('C1',C1,'C2',C2,'C3',C3), ...
                       'fileCSV', csvPath);
            axInfo.scatter = [axInfo.scatter; s]; %#ok<AGROW>
        end

        extracted.axes = [extracted.axes; axInfo]; %#ok<AGROW>
    end

    % Save everything to MAT
    matPath = fullfile(outDir, sprintf('%s_extracted_%s.mat', base, timestamp));
    save(matPath, 'extracted');

    % Build summary
    summary = struct();
    summary.figFile     = extracted.sourceFig;
    summary.outputDir   = outDir;
    summary.savedCSV    = savedFiles(:);
    summary.matFile     = matPath;
    summary.numAxes     = numel(extracted.axes);
    summary.numLines    = sum(arrayfun(@(a) numel(a.lines), extracted.axes));
    summary.numScatters = sum(arrayfun(@(a) numel(a.scatter), extracted.axes));

    fprintf('Extracted from: %s\n', summary.figFile);
    fprintf('Axes: %d | Lines: %d | Scatters: %d\n', summary.numAxes, summary.numLines, summary.numScatters);
    fprintf('CSV files saved to: %s\n', summary.outputDir);
    fprintf('MAT saved to: %s\n', summary.matFile);
end

% ===== helpers =====
function v = tryGet(h, prop, defaultVal)
    if isgraphics(h) && isprop(h, prop)
        v = get(h, prop);
    else
        v = defaultVal;
    end
end

function t = tryAxTitle(ax)
    t = "";
    try
        T = get(ax,'Title');
        if isgraphics(T) && isprop(T,'String'), t = string(T.String); end
    catch
    end
end

function tf = isLegendOn(h)
    tf = true;
    try
        L = h.Annotation.LegendInformation;
        if isprop(L,'IconDisplayStyle')
            tf = strcmpi(L.IconDisplayStyle,'on');
        end
    catch
        % default true
    end
end

function tf = isErrorbarChild(h)
    tf = false;
    try
        hg = ancestor(h, 'hggroup');
        if ~isempty(hg)
            tg = get(hg,'Tag');
            if ischar(tg), tf = contains(lower(tg),'errorbar'); end
        end
    catch
    end
end

function tf = hasAtLeastTwoFinitePoints(X, Y)
    X = X(:); Y = Y(:);
    m = min(numel(X), numel(Y));
    if m < 2, tf = false; return; end
    ok = isfinite(X(1:m)) & isfinite(Y(1:m));
    tf = nnz(ok) >= 2;
end

function objs = dedupByData(objs, kind)
    if isempty(objs), return; end
    dup = false(size(objs));
    for a = 1:numel(objs)
        if dup(a), continue; end
        [Xa,Ya,Za] = getXYZ(objs(a));
        for b = a+1:numel(objs)
            if dup(b), continue; end
            [Xb,Yb,Zb] = getXYZ(objs(b));
            sameXY = isequaln(Xa, Xb) && isequaln(Ya, Yb);
            if sameXY
                if ~isempty(Za) || ~isempty(Zb)
                    if isequaln(Za, Zb)
                        dup(b) = true;
                    end
                else
                    dup(b) = true;
                end
            end
        end
    end
    objs = objs(~dup);
end

function [X, Y, Z] = getXYZ(h)
    X = tryGet(h,'XData',[]);
    Y = tryGet(h,'YData',[]);
    Z = [];
    if isprop(h,'ZData')
        Z = tryGet(h,'ZData',[]);
        if isempty(Z) || all(~isfinite(Z)), Z = []; end
    end
    % ensure vectors
    X = X(:); Y = Y(:);
    if ~isempty(Z), Z = Z(:); end
end

function T = makeLineTable(X, Y, Z)
    T = table();
    if ~isempty(X), T.X = X; end
    if ~isempty(Y), T.Y = Y; end
    if ~isempty(Z) && numel(Z) == height(T), T.Z = Z; end
end

function [C1,C2,C3] = normalizeC(C, n)
    C1 = []; C2 = []; C3 = [];
    if isempty(C), return; end
    if isscalar(C)
        C1 = repmat(C, n, 1);
        return;
    end
    if ismatrix(C) && size(C,1) == n && size(C,2) == 3
        C1 = C(:,1); C2 = C(:,2); C3 = C(:,3);
        return;
    end
    C = C(:);
    if numel(C) < n, C(end+1:n,1) = NaN; end
    if numel(C) > n, C = C(1:n); end
    C1 = C;
end

function S = normalizeSize(S, n)
    if isempty(S), return; end
    if isscalar(S), S = repmat(S, n, 1); return; end
    S = S(:);
    if numel(S) < n, S(end+1:n,1) = NaN; end
    if numel(S) > n, S = S(1:n); end
end

function safeClose(fig)
    if ishghandle(fig)
        close(fig);
    end
end
