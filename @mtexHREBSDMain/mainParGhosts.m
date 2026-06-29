function mainParGhosts(obj, input, options)
    pool = gcp("nocreate");
    if isempty(pool)
        pool = parpool(options.numCores);
    end
    wb = waitbar(0, "Starting analysis...");
    cleanupObj = onCleanup(@()close(wb));
    ebsd = obj.ebsdSubset(input, options.Subset);
    subsetIds = ebsd.id;
    totalPoints = length(subsetIds);
    numChunks = totalPoints/options.ChunkSize;
    cfs = options.calcFStrain;
    cfa = options.calcFAlpha;
    [skipder, skipstrain] = initializeSkip(input, obj.doStrain, obj.doAlpha);
    [out, start, kStart] = initializeOutput(totalPoints, options);
    waitbar(0, wb, "Beginning analysis...");
    nUpdateWaitbar(numChunks, wb, start);
    k = kStart;
    tic
    for i = start:numChunks
        skipe = zeros(1,options.ChunkSize);
        skipd = zeros(1,options.ChunkSize);
        for l = 1:options.ChunkSize
            indexIter = (i-1)*options.ChunkSize + l;
            iditer = subsetIds(indexIter);
            skipe(l) = skipstrain(iditer);
            skipd(l) = skipder(iditer);
            ids = patPackageIds(input, iditer);
            try
                pats(l,:) = patPackage(input, ids, skipe(l), skipd(l));
            catch me
                obj.errorLog{end+1} = {indexIter, me};
            end
        end
        % for troubleshooting Main()
%         outIter = Main(cfs, cfa, pats, skipd, obj.doStrain, obj.doAlpha) 
        if ~exist('F', 'var')
            F(1) = parfeval(pool, @Main, 1,...
                cfs, cfa, pats, skipe, skipd...
                );
        else
            F(end+1) = parfeval(pool, @Main, 1,...
                cfs, cfa, pats, skipe, skipd...
                );
        end
        if length(F) >= options.QueueSizeMult*options.numCores
            if options.Verbose
                clc
                disp("Chunk # "+num2str(i));
                disp(F)
            end
            wait(F(1:options.numCores))
        end
        if sum(strcmp({F.State}, 'finished')) %&& length(F) > 2*options.numCores
            if options.Verbose
                clc
                disp("Chunk # "+num2str(i));
                disp(F)
            end
            for j = 1:sum(strcmp({F.State}, 'finished'))
                [idx, outIter, k, range, errormsg] = readFutures( ...
                    F, k, options.ChunkSize ...
                    );
                if ~isempty(errormsg)
                    obj.errorLog{end+1} = {k, errormsg};
                end
                out(range, 1:49) = outIter;
                out(range, 50) = k;
                F(idx) = [];
            end    
        end
        if ~isempty(options.Autosave)
            if ~mod(i,options.AutosaveInterval)
                disp("Autosaving to "+options.Autosave+"...")
                save(options.Autosave, "out")
            end
        end
    end
    if options.Verbose
        clc
        disp(F)
    end
    wait(F)
    for j = 1:length(F)
        [idx, outIter, k, range, errormsg] = readFutures( ...
            F, k, options.ChunkSize ...
            );
        if ~isempty(errormsg)
            obj.errorLog{end+1} = {k, errormsg};
        end
        out(range, 1:49) = outIter;
        out(range, 50) = k;
        F(idx) = [];
        nUpdateWaitbar;
    end
    if ~isempty(options.Autosave)
        if ~mod(i,options.AutosaveInterval)
            disp("Autosaving to "+options.Autosave+"...")
            save(options.Autosave, "out")
        end
    end
    obj.readOutput(out);
%     close(wb)
    toc
end


function scanIndicies = getScanIndicies(F, chunkSize)
    scanIndicies = zeros(1,chunkSize);
    for i = 1:chunkSize
        pat = F.InputArguments{1,3}(i,1);
        scanIndicies(i) = pat.scanIndex;
    end
end



function [out, start, kStart] = initializeOutput(totalPoints, options)
    if ~isempty(options.Autosave) && isfile(options.Autosave)
        out = load(options.Autosave).out;
        if ~isempty(options.Subset)
            maxGlobalId = max(out(:,1));
            subsetId = find(options.Subset == maxGlobalId);
            start = subsetId/options.ChunkSize;
        else
            start = max(out(:,1))/options.ChunkSize;
        end
        kStart = max(out(:,end));
    else 
        out = zeros(totalPoints, 450);
        start = 1;
        kStart = 0;
    end
end


function [i, out, k, range, errorMsg] = readFutures(F, k, chunk)
    % check for errors
    try
        [i, out] = fetchNext(F);
        errorMsg = [];
    catch me
        if strcmp(me.identifier, ...
                'MATLAB:parallel:future:FetchNextFutureErrored' ...
                )
            errorMsg = me;
            hasError = ~cellfun(@isempty, {F.Error});
            i = find(hasError, 1, 'first');
            out = zeros(chunk,49);
            for m = 1:chunk
                pat = F(i).InputArguments{1,3}(m,1);
                out(1,m) = pat.scanIndex;
            end
        end
    end
    k = k + 1;
    lo = (k-1)*chunk + 1;
    hi = (k-1)*chunk + size(out,1);
    range = lo:hi;
    nUpdateWaitbar;
end
    


function [skipDerivs, skipStrain] = initializeSkip(input, doStrain, doAlpha)
    nx = input.scan.Nx;
    ny = input.scan.Ny;
    if doAlpha
        xind = round(input.ebsd.prop.x/input.scan.xStep + 1);
        yind = round(input.ebsd.prop.y/input.scan.yStep + 1);
        skipdx = xind == nx;
        skipdy = yind == ny;
        skipDerivs = skipdx + skipdy;
    else
        skipDerivs = ones(1,length(input.ebsd));
    end
    if doStrain
        skipStrain = zeros(1,length(input.ebsd));
    else
        skipStrain = ones(1,length(input.ebsd));
    end

end


function out = Main(cfs, cfa, pats, skipStrain, skipDer)
    if size(pats,1) > 1
        out = zeros(size(pats,1), 46);
        for i = 1:size(pats,1)
            out(i,:) = Main( ...
                cfs, cfa, pats(i,:), skipStrain(i), skipDer(i) ...
                );
        end
    else
        out = zeros(1,49);
        out(1) = pats(1).scanIndex;
        if ~skipStrain
            cfs.getBeta(pats(1), pats(2));
            out(2:10) = cfs.beta;
            out(11:19) = cfs.g;
            out(38:40) = pats(1).patternCenter;
            if cfs.AnalysisType == 2
                out(47:49) = cfs.patternCenter;
            end
            out(41) = cfs.SSE;
            out(42) = cfs.R2;
        end 
        if ~skipDer
            cfa.getBeta(pats(3), pats(1));
            out(20:28) = cfa.beta;
            out(43) = cfa.SSE;
            out(44) = cfa.R2;
            cfa.getBeta(pats(4), pats(1));
            out(29:37) = cfa.beta;
            out(45) = cfa.SSE;
            out(46) = cfa.R2;
        end
    end
end


function pats = patPackage(input, ids, skipStrain, skipDer)
    pats(1) = input.get_pattern(ids(1));
    if ~skipStrain
        pats(2) = input.get_pattern(ids(2));
    else
        pats(2) = pats(1); % Fill with test pattern to avoid reading pat
    end
    if ~skipDer
        pats(3) = input.get_pattern(ids(3));
        pats(4) = input.get_pattern(ids(4));
    else
        pats(3) = pats(1);
        pats(4) = pats(1);
    end
end


function ids = patPackageIds(input, i)
    ids = zeros(1,4);
    ids(1) = i;
    ids(2) = input.refIds(input.ebsd.grainId(i));
    [ids(3), ids(4)] = getDerivIds(input,i);
end


function [idx,idy] = getDerivIds(input, i)
    [~,~,ext] = fileparts(input.scan.scanfile);
    if any(strcmp(ext, [".ctf", ".h5oina"]))
        nx = input.scan.Nx;
        ny = input.scan.Ny;
        [x,y] = ind2sub([nx,ny], i);
        if x+1 > nx || y+1 > ny
            idx = 0;
            idy = 0;
        else
            idx = sub2ind([nx,ny], x+1, y);
            idy = sub2ind([nx,ny], x, y+1);
        end
    else % Double check this weirdness isn't happening with oxford systems...
        nx = input.scan.Ny;
        ny = input.scan.Nx;
        [x,y] = ind2sub([nx,ny], i);
        if x+1 > nx || y+1 > ny
            idx = 0;
            idy = 0;
        else
            idx = sub2ind([nx,ny], x, y+1);
            idy = sub2ind([nx,ny], x+1, y);
        end
    end
end