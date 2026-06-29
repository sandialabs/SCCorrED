function obj = run_grainClassic(obj, mtexHREBSD)
    disp('Classic')
    tic
    numPoints = length(mtexHREBSD.ebsd);
%     obj.globalIds = zeros(numPoints,1);
    localRefId = find(mtexHREBSD.ebsd.id == mtexHREBSD.refIds);
    refPat = mtexHREBSD.get_pattern(localRefId);
    p = waitbar(0/numPoints, 'Starting classic analysis');
    for j = 1:numPoints
        try
            testPat = mtexHREBSD.get_pattern(j);
%             obj.globalIds(j) = testPat.scanIndex;
            calcF = classicCalcF(refPat, testPat, mtexHREBSD);
            obj.beta(:,:,j) = calcF.beta;
            obj.F(:,:,j) = calcF.F;
            obj.g(:,:,j) = calcF.g;
            obj.fitMetrics.SSE(calcF.j) = calcF.fit.metrics.SSE;
            obj.fitMetrics.R2(calcF.j) = calcF.fit.metrics.R2;
        catch me
        end
        if ~mod(j, 10)
            time_j = toc;
            obj.update_waitbar(time_j, numPoints, j, p)
        end
    end
    obj.beta = obj.beta(:,:,1:numPoints);
    obj.F = obj.F(:,:,1:numPoints);
    obj.g = obj.g(:,:,1:numPoints);

%     obj.fitMetrics.SSE(calcF.j) = calcF.fit.metrics.SSE;
%     obj.fitMetrics.R2(calcF.j) = calcF.fit.metrics.R2;
    time_j = toc;
    obj.update_waitbar(time_j, numPoints, j, p)
end