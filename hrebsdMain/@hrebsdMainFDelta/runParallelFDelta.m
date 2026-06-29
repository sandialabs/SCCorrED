function obj = runParallelFDelta(obj, mtexHREBSD)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ebsd = obj.ebsdSubset(mtexHREBSD, obj.args);
    wb = waitbar(0, "Starting analysis... | Slicing data");
    totalPoints = length(ebsd);
    grainsToRun = unique(ebsd.grainId);
    inputs = mtexHREBSD_lite(mtexHREBSD);
    [reorderIds, ~] = getIdOrder(ebsd, grainsToRun);
    waitbar(0, wb, "Setting up processing...");
    F = zeros(3,3,totalPoints);
    g = zeros(3,3,totalPoints);
    res = zeros(3, totalPoints);
    pcs = zeros(3, totalPoints);
%     SSE = zeros(totalPoints, 1);
%     R2 = zeros(totalPoints, 1);
    ph = zeros(totalPoints, 1);
    logMessages = cell(totalPoints, 1);
%     nUpdateWaitbar(totalPoints, wb);
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    waitbar(0/2, wb, "Sending input object to workers...");
    C = parallel.pool.Constant(mtexHREBSD);
    waitbar(2/2, wb, "Beginning analysis...");
    nUpdateWaitbar(totalPoints, wb);
    tic
    parfor i = 1:length(reorderIds)
        testId = reorderIds(i);
%         refId = reorderRefIds(i);
        testPat = C.Value.get_pattern(testId);
        try
%             [newPat, newF, resi] = FDelta(testPat, C.Value);
            calcF = FDeltaCalcF(testPat, mtexHREBSD, varargin)
%             calcF = classicCalcF(refPatsConst.Value{refId}, testPat, inputs);
            F(:,:,i) = calcF.F;
            g(:,:,i) = calcF.g;
%             SSE(i) = calcF.fit.metrics.SSE;
%             R2(i) = calcF.fit.metrics.SSE;
%             res(:,i) = resi;
            pcs(:,i) = calcF.patternCenter;
            ph(i) = calcF.peakHeight;
        catch me
            me
            logMessages{i} = {i, me};
        end
        send(D,1)
    end
    obj.F(:,:,reorderIds) = F;
    obj.beta = obj.F - eye(3);
    obj.g(:,:,reorderIds) = g;
    obj.residual(:,reorderIds) = res;
    obj.pcs(:,reorderIds) = pcs;
    obj.peakHeights(reorderIds) = ph;
%     obj.fitMetrics.SSE(reorderIds) = SSE;
%     obj.fitMetrics.R2(reorderIds) = R2;
    for i = 1:totalPoints
        msgIter = logMessages{i};
        if ~isempty(msgIter)
            obj.logMessage{end+1} = msgIter;
        end
    end
%     obj.logMessage{reorderIds} = logMessages;
%     close(wb)
    toc
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
end