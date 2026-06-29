function obj = run_grainHybrid(obj, mtexHREBSD)
    tic
    numPoints = length(mtexHREBSD.ebsd);
    obj.globalIds = zeros(numPoints,1);
    localRefId = find(mtexHREBSD.ebsd.id == mtexHREBSD.refIds);
    refPat = mtexHREBSD.get_pattern(localRefId);
    refF = dynamicSimulatedCalcF(refPat, mtexHREBSD, 'ConvergencePlots',1);
    p = waitbar(0/numPoints, 'Starting hybrid analysis');
    for j = 1:numPoints
        testPat = mtexHREBSD.get_pattern(j);
        obj.globalIds(j) = testPat.scanIndex;
        calcF = classicCalcF(refPat, testPat, mtexHREBSD);
        diffF = calcF.F + refF.F;
        obj.beta(:,:,j) = calcF.beta;
        obj.F(:,:,j) = calcF.F - diffF;
        if ~mod(j, 100)
            time_j = toc;
            obj.update_waitbar(time_j, numPoints, j, p)
        end
    end
    time_j = toc;
    obj.update_waitbar(time_j, numPoints, j, p)
end