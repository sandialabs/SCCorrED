function obj = run_analysis(obj, mtexHREBSD)
%run_analysis Summary of this function goes here
%   Detailed explanation goes here

    disp("Starting single process HREBSD analysis...")
    tic
    grains = unique(mtexHREBSD.ebsd.grainId);
    numGrains = length(grains);
%     numGrains = length(mtexHREBSD.grains);
%     numGrains = 25;
    waitStr = "Starting single process HREBSD analysis for " ...
                + num2str(numGrains) + " grains";
    waitbar(0/length(numGrains), obj.progressBar, waitStr);
    for i = 1:numGrains % looping through grains
        grainToUse = grains(i);
        grain_i = mtexHREBSD.get_grain(grainToUse);
        if isnan(grain_i.ebsd.rotations.phi1)
            continue
        end
        numPoints = length(grain_i.ebsd);
%         viableRefs = viableRefIds(grain_i.ebsd, mtexHREBSD.analysis.numRefIds);
        refIdGrain = find(grain_i.ebsd.id == mtexHREBSD.refIds(grainToUse),1);
        refPattern = grain_i.get_pattern(refIdGrain);
        for j = 1:numPoints % looping though points in grain
            testPattern = grain_i.get_pattern(j);
            try
                calcF_j = feval(obj.methodCalcF, refPattern, testPattern, grain_i);
                obj.beta(:,:,testPattern.scanIndex) = calcF_j.beta;
                obj.F(:,:,testPattern.scanIndex) = calcF_j.F;
            catch ME
                disp("Failed at " + num2str(i) + "th grain at the " +num2str(j)+ "th point...")
                disp(ME.message)
            end
        end
        time_i = toc;
        [h,m,s] = hms(seconds(time_i));
        waitStr = "Last completed grain " + num2str(i) + ...
            ", elapsed time " + num2str(h) + ":" + num2str(m) + ...
            ":" + num2str(s);
        waitbar(i/numGrains, obj.progressBar, waitStr);
    end
%     toc
    obj.completionTime = toc;
end

function check = check_samePattern(refPattern, testPattern)
    if refPattern.scanIndex == testPattern.scanIndex
        check = 1;
    else
        check = 0;
    end
end
