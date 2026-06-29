function obj = runParallel3(obj, mtexHREBSD)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    wb = waitbar(0, "Starting analysis... ");
    totalPoints = length(mtexHREBSD.ebsd);
    grainsToRun = mtexHREBSD.grains.id;
%     refPats = getRefPats(mtexHREBSD, grainsToRun, wb);
    waitbar(0, wb, "Setting up processing...");
    beta = zeros(3,3,totalPoints);
    g = zeros(3,3,totalPoints);
    SSE = zeros(totalPoints, 1);
    R2 = zeros(totalPoints, 1);
%     logMessages = cell(totalPoints, 1);
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    waitbar(0/1, wb, "Sending input object to workers...");
    C = parallel.pool.Constant(mtexHREBSD);
    waitbar(1/1, wb, "Sending input object to workers...");
    nUpdateWaitbar(totalPoints, wb);
    for i = 1:length(grainsToRun)
        inds = find(mtexHREBSD.ebsd.grainId == grainsToRun(i));
        if length(inds) < 2
            continue
        else
            refPatConst = parallel.pool.Constant( ...
                mtexHREBSD.get_pattern(mtexHREBSD.refIds(grainsToRun(i))) ...
                );
            betaTemp = zeros(3,3,length(inds));
            gTemp = zeros(3,3,length(inds));
            SSETemp = zeros(length(inds), 1);
            R2Temp = zeros(length(inds), 1);
            if i == 1
                tic
            end
            parfor j = 1:length(inds)
                testInd = inds(j);
                testPat = C.Value.get_pattern(testInd);
                calcF = classicCalcF(refPatConst.Value, testPat, C.Value);
                betaTemp(:,:,j) = calcF.beta;
                gTemp(:,:,j) = calcF.g;
                SSETemp(j) = calcF.fit.metrics.SSE;
                R2Temp(j) = calcF.fit.metrics.R2;
            end
            for k = 1:length(inds)
                beta(:,:,inds) = betaTemp;
                g(:,:,inds) = gTemp;
                SSE(inds) = SSETemp;
                R2(inds) = R2Temp;
            end
            send(D,1)
        end  
    end
end


% function refPats = getRefPats(mtexHREBSD, grainsToRun, wb)
%     refPats = cell(length(grainsToRun), 1);
%     for i = 1:length(grainsToRun)
%         waitbar(i/length(grainsToRun), wb, "Getting reference patterns...");
%         refPats{i} = mtexHREBSD.get_pattern(mtexHREBSD.refIds(grainsToRun(i)));
%     end
% end
