function mainSingle(obj, input, options)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% Updated to include R2, SSE, etc..., TJB, 2024-11-22
    ebsd = obj.ebsdSubset(input, options.Subset);
    cfs = options.calcFStrain;
    cfa = options.calcFAlpha;
    wb = waitbar(0, "Starting analysis... | Slicing data");
    cleanupObj = onCleanup(@()close(wb));
    grainsToRun = input.grains.id;
    [reorderIds, reorderRefIds] = getIdOrder(ebsd, grainsToRun);
    totalPoints = length(reorderIds);
    logger = {};
    beta = zeros(3,3,totalPoints);
    dbx = zeros(3,3,totalPoints);
    dby = zeros(3,3,totalPoints);
    g = zeros(3,3,totalPoints);
    R2 = zeros(1,totalPoints);
    SSE = zeros(1,totalPoints);
    pcinitial = zeros(3, totalPoints);
    if cfs.AnalysisType == 2
        pccalibrated = zeros(3,totalPoints);
    end
    waitbar(0, wb, "Beginning analysis...");
    nUpdateWaitbar(totalPoints, wb);
    tic
    for i = 1:length(reorderIds)
        testId = reorderIds(i);
        pt = input.get_pattern(testId);
        if obj.doStrain
            try
                strainAnalysis(cfs, input, pt, reorderRefIds, i);
                beta(:,:,i) = cfs.beta;
                g(:,:,i) = cfs.g;
                R2(i) = cfs.R2;
                SSE(i) = cfs.SSE;
                pcinitial(:,i) = pt.patternCenter;
                if cfs.AnalysisType == 2
                    pccalibrated(:,i) = cfs.patternCenter;
                end
            catch me
                logger{end+1} = {testId, me};
            end
        end
        if obj.doAlpha
            try
                [dbx(:,:,i),dby(:,:,i)] = alphaAnalysis(cfa, input, pt, i);
            catch me
                logger{end+1} = {testId, me};
            end
        end
        nUpdateWaitbar;
    end
    if obj.doStrain
        obj.beta(:,:,reorderIds) = beta;
        obj.F = obj.beta + eye(3);
        obj.g(:,:,reorderIds) = g;
        obj.R2(reorderIds) = R2;
        obj.SSE(reorderIds) = SSE;
        obj.PCInitial(:,reorderIds) = pcinitial;
        if cfs.AnalysisType == 2
            obj.PCCalibrated(:,reorderIds) = pccalibrated;
        end
    end
    if obj.doAlpha
        obj.dbetadx(:,:,reorderIds) = dbx;
        obj.dbetady(:,:,reorderIds) = dby;
    end
end


function [dx,dy] = alphaAnalysis(cfa, input, pt, i)
    nx = input.scan.Nx;
    ny = input.scan.Ny;
    xstep = input.scan.xStep;
    ystep = input.scan.yStep;
    [x,y] = ind2sub([nx,ny], i);
    if x ~= nx && y ~= ny
        idx = sub2ind([nx,ny], x+1, y);
        idy = sub2ind([nx,ny], x, y+1);
        px = input.get_pattern(idx);
        py = input.get_pattern(idy);
        cfa.getBeta(px, pt);
        dx = cfa.beta/(xstep*1e-6);
        cfa.getBeta(py,pt);
        dy = cfa.beta/(ystep*1e-6);
    else
        dx = zeros(3,3);
        dy = zeros(3,3);
    end
end


function strainAnalysis(cfs, input, pt, reorderRefIds, i)
    if ~cfs.AnalysisType
        refId = input.refIds(reorderRefIds(i));
        pr = input.get_pattern(refId);
        cfs.getBeta(pt,pr);
    else 
        cfs.getBeta(pt);
    end
end


function [reorderIds, reorderRefIds] = getIdOrder(ebsd, subsetGrains)
    reorderIds = zeros(length(ebsd), 1);
    reorderRefIds = zeros(length(ebsd), 1);
    currIndex = 1;
    for i = 1:length(subsetGrains)
        idsInGrain = ebsd.id(ebsd.grainId == subsetGrains(i));
        numIds = length(idsInGrain);
        reorderIds(currIndex:currIndex+numIds-1) = idsInGrain;
        reorderRefIds(currIndex:currIndex+numIds-1) = i;
        currIndex = currIndex + numIds;
    end
    reorderIds(reorderIds==0) = [];
    reorderRefIds(reorderIds==0) = [];
end