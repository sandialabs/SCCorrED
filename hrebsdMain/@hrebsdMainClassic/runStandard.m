function obj = runStandard(obj, mtexHREBSD)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ebsd = obj.ebsdSubset(mtexHREBSD, obj.args);
    wb = waitbar(0, "Starting analysis... | Slicing data");
%     totalPoints = length(ebsd);
%     grainsToRun = unique(mtexHREBSD)
    grainsToRun = mtexHREBSD.grains.id;
    [reorderIds, reorderRefIds] = getIdOrder(ebsd, grainsToRun);
    totalPoints = length(reorderIds);
    refPats = cell(length(grainsToRun),1);
    for i = 1:length(grainsToRun)
        waitbar(i/length(grainsToRun), wb, "Getting reference patterns...");
        refPats{i} = mtexHREBSD.get_pattern(mtexHREBSD.refIds(grainsToRun(i)));
    end
    waitbar(0, wb, "Setting up processing...");
    beta = zeros(3,3,totalPoints);
    g = zeros(3,3,totalPoints);
    SSE = zeros(totalPoints, 1);
    R2 = zeros(totalPoints, 1);
    ph = zeros(totalPoints, 1);
    logMessages = cell(totalPoints, 1);
%     nUpdateWaitbar(totalPoints, wb);
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    waitbar(0/2, wb, "Sending input object to workers...");
%     C = parallel.pool.Constant(mtexHREBSD);
    waitbar(1/2, wb, "Sending reference patterns to workers...");
%     refPatsConst = parallel.pool.Constant(refPats);
    waitbar(2/2, wb, "Beginning analysis...");
    nUpdateWaitbar(totalPoints, wb);
    tic
    for i = 1:length(reorderIds)
        testId = reorderIds(i);
        refId = reorderRefIds(i);
        testPat = mtexHREBSD.get_pattern(testId);
%         testPat = mtexHREBSD.get_pattern(testId);
        try
%             disp(refPatsConst.Value{refId})
            calcF = classicCalcF(refPats{refId}, testPat, mtexHREBSD);
            beta(:,:,i) = calcF.beta;
            g(:,:,i) = calcF.g;
            SSE(i) = calcF.fit.metrics.SSE;
            R2(i) = calcF.fit.metrics.R2;
            ph(i) = calcF.peakHeight;
        catch me
            disp(me)
            logMessages{i} = {i, me};
        end
        send(D,1)
    end
    obj.beta(:,:,reorderIds) = beta;
    obj.F = obj.beta + eye(3);
    obj.g(:,:,reorderIds) = g;
    obj.fitMetrics.SSE(reorderIds) = SSE;
    obj.fitMetrics.R2(reorderIds) = R2;
    obj.peakHeights(reorderIds) = ph;
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
    reorderIds(reorderIds==0) = [];
    reorderRefIds(reorderIds==0) = [];
end