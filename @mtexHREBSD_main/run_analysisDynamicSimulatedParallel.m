function obj = run_analysisDynamicSimulatedParallel(obj, mtexHREBSD)
    tic
    numPoints = length(mtexHREBSD.ebsd);
%     numPoints = 8;
    pool = gcp("nocreate");
    if isempty(pool)
        parpool(mtexHREBSD.analysis.numCores);
        pool = gcp("nocreate");
    end
    N = pool.NumWorkers;
    h = waitbar(0, "Settings up parallel job...");
    timePerPoint = 0;
    for j = 1:N:numPoints
        for i = 0:N-1
            if j + i <= numPoints
                temp = mtexHREBSD;
                testPat = mtexHREBSD.get_pattern(j+i);
                if j+i == 1
                    futures(1) = parfeval(@dynamicSimulatedCalcF, 1, testPat, temp);
                else
                    futures(i+1) = parfeval(@dynamicSimulatedCalcF, 1, testPat, temp);
                end
%                 if mod(length(futures), N) == 0
%                     wait(futures)
%                 end
                if i == N-1
                    wait(futures)
                end
            else
                break
            end
        end
        for i = 0:N-1
            if j + i <= numPoints
                try
                    [~,calcF] = fetchNext(futures);
                    ind = calcF.scanIndex;
                    obj.F(:,:,ind) = calcF.F;
                    obj.beta(:,:,ind) = calcF.beta;
                    obj.g(:,:,ind) = calcF.g;
                    obj.fitMetrics.SSE(ind) = calcF.fit.metrics.SSE;
                    obj.fitMetrics.R2(ind) = calcF.fit.metrics.R2;
                catch ME
                    obj.logMessage{end+1} = {j+i, ME};
                end
                runTime = toc;
                timePerPoint = (timePerPoint + runTime)/(j+i);
                message = get_progressMessage((j+i), runTime, timePerPoint, numPoints);
                waitbar((j+i)/numPoints, h, message)
            else
                break
            end
        end
        obj.futures = futures;
    end
    obj.completionTime = toc;
end


function message = get_progressMessage(i, runTime, timePerPoint, numPoints)
    [h,m,s] = hms(seconds(runTime));
    elapsedTime = "Elapsed time: "+num2str(h)+":"+num2str(m)+":"+num2str(s);
    estimatedTimeRemaining = (numPoints - i)*timePerPoint;
    [hr,mr,sr] = hms(seconds(estimatedTimeRemaining));
    estimate = "Estimate time remaining: "+num2str(hr)+":"+num2str(mr)+":"+num2str(sr);
    message = ["Points completed: " + num2str(i) + "/" + num2str(numPoints) + ", " + elapsedTime, estimate];
end