function obj = run_analysisHybridParallel(obj, mtexHREBSD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    numGrains = length(mtexHREBSD.grains);
%     numGrains = 20;
    nUpdateWaitbar(numGrains, obj.progressBar)
    storageCell = cell(numGrains,1);
    refF = cell(numGrains, 1);
    pool = gcp("nocreate");
    if isempty(pool)
        parpool(mtexHREBSD.analysis.numCores);
    end
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar);
    tic
    parfor i = 1:numGrains
        temp = mtexHREBSD;
        grain_i = temp.get_grain(i);
        numPoints = length(grain_i.ebsd);
        calcFs = cell(1,numPoints);
        refIdGrain = find(grain_i.ebsd.id == temp.refIds(i));
        refPattern = grain_i.get_pattern(refIdGrain);
        if ~isnan(refPattern.rotations(1))
            refF{i} =  dynamicSimulatedCalcF(refPattern, grain_i);
            for j = 1:numPoints
                if ~isnan(grain_i.phi1(j))
                    testPattern = grain_i.get_pattern(j);
                    try
                        calcF_ij = classicCalcF(refPattern, testPattern, grain_i);
                        % subtract the difference between iteration
                        % calcF_ij and the strain in the reference pattern
                        % to get the absolute strain!
                        calcF_ij.F = calcF_ij.F*refF{i}.F;
                        calcF_ij.beta = calcF_ij.F - eye(3);
                        calcFs{j} = calcF_ij;
                    catch ME
                        disp("Failed at " + num2str(i) + "th grain at the " +num2str(j)+ "th point...")
                        disp(ME.message)
                    end
                end
            end
            storageCell{i,1} = calcFs;
        end
        send(D,1);
    end
    obj.refF = refF;
%     obj.out = storageCell;
    [obj.F, obj.beta, obj.g, obj.fitMetrics] = get_dataFromStorageCell(storageCell, mtexHREBSD.scan.scanLength);
    obj.completionTime = toc;
end


function [F, beta, g, fitMetrics] = get_dataFromStorageCell(storageCell, scanLength)
    F = zeros(3,3,scanLength);
    beta = zeros(3,3,scanLength);
    g = zeros(3,3,scanLength);
    fitMetrics = struct("SSE", zeros(1,scanLength), "R2", zeros(1,scanLength));
    for i = 1:length(storageCell)
        entries = storageCell{i,1};
        for j = 1:length(entries)
            if ~isempty(entries{j})
                index = entries{j}.scanIndex;
                F(:,:,index) = entries{j}.F;
                beta(:,:,index) = entries{j}.beta;
                g(:,:,index) = entries{j}.g;
                fitMetrics.SSE(index) = entries{j}.fit.metrics.SSE;
                fitMetrics.R2(index) = entries{j}.fit.metrics.R2;
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