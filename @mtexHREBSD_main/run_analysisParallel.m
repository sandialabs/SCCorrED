function obj = run_analysisParallel(obj,mtexHREBSD)
%run_analysisParallel Summary of this function goes here
%   Detailed explanation goes here

    numGrains = length(mtexHREBSD.grains);
%     numGrains = 15;
    nUpdateWaitbar(numGrains, obj.progressBar)
%     methodCalcF = obj.methodCalcF;
    numWorkers = mtexHREBSD.analysis.numCores;
    storageCell = cell(numGrains,1);
    pool = gcp("nocreate");
    if isempty(pool)
        parpool(mtexHREBSD.analysis.numCores);
    end
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar);
    tic
    parfor i = 1:numGrains % looping through grains
        temp = mtexHREBSD;
        grain_i = temp.get_grain(i);
%         if ~isnan(grain_i.phi1)
        numPoints = length(grain_i.ebsd);
%         testIds = zeros(1,numPoints);
        calcFs = cell(1,numPoints);
%     %         viableRefs = viableRefIds(grain_i.ebsd, mtexHREBSD.analysis.numRefIds);
        refIdGrain = find(grain_i.ebsd.id == temp.refIds(i));
        refPattern = grain_i.get_pattern(refIdGrain);
        for j = 1:numPoints % looping though points in grain
            if ~isnan(grain_i.phi1(j))
                testPattern = grain_i.get_pattern(j);
%                 testInd = testPattern.scanIndex;
%                 testIds(j) = testInd;
                try
                    calcFs{j} = classicCalcF(refPattern, testPattern, grain_i);
                catch ME
                    disp("Failed at " + num2str(i) + "th grain at the " +num2str(j)+ "th point...")
                    disp(ME.message)
                end
            end
        end
        storageCell{i,1} = calcFs;
        send(D, 1);
    end
%     obj.out = storageCell;
    [obj.F, obj.beta, obj.g] = get_dataFromStorageCell(storageCell, mtexHREBSD.scan.scanLength);
    obj.completionTime = toc;
end


function [F, beta, g] = get_dataFromStorageCell(storageCell, scanLength)
    F = zeros(3,3,scanLength);
    beta = zeros(3,3,scanLength);
    for i = 1:length(storageCell)
        entries = storageCell{i,1};
        for j = 1:length(entries)
            if ~isempty(entries{j})
                index = entries{j}.scanIndex;
                F(:,:,index) = entries{j}.F;
                beta(:,:,index) = entries{j}.beta;
                g(:,:,index) = entries{j}.g;
            end
        end
    end
end
    


function check = check_samePattern(refPattern, testPattern)
    if refPattern.scanIndex == testPattern.scanIndex
        check = 1;
    else
        check = 0;
    end
end


function p = nUpdateWaitbar(data, h)
    persistent TOTAL COUNT H
    if nargin == 2
        H = h;
        TOTAL = data;
        COUNT = 0;
    else
        COUNT = 1 + COUNT;  
        p = COUNT/TOTAL;
        t = toc;
        [h,m,s] = hms(seconds(t));
        waitStr = "Completed grains " + num2str(COUNT) + ...
            ", elapsed time " + num2str(h) + ":" + num2str(m) + ...
            ":" + num2str(s);
        waitbar(p, H, waitStr)
    end
end