classdef mtexHREBSD_main
    %mtexHREBSD_main Summary of this class goes here
    %   Detailed explanation goes here

    properties
        F
        beta
        g
        fitMetrics
        alpha
        completionTime
        version
        refF
        logMessage
        futures
    end

    properties(Hidden=true)
%         globalIds
%         out
        methodCalcF
        methodAnalysis
        progressBar
        ft % = frameTransformations
    end

    methods
        function obj = mtexHREBSD_main(mtexHREBSD, varargin)
            %mtexHREBSD_main Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj.initialize(mtexHREBSD);
            obj.methodAnalysis = obj.get_methodAnalysis(mtexHREBSD);
%             obj.methodCalcF = obj.get_methodCalcF(mtexHREBSD);
            if nargin == 1
%                 obj.methodAnalysis = obj.get_methodAnalysis(mtexHREBSD);
                if mtexHREBSD.isGrain
                    obj = obj.run_grain(mtexHREBSD);
                else
                    if mtexHREBSD.analysis.numCores > 1
                        obj = obj.run_parallel(mtexHREBSD);
                    else 
                        obj = obj.run_notParallel(mtexHREBSD);
                    end
                end
            end
        end


        function E = get_strain(obj)
            E = 0.5*(obj.beta + permute(obj.beta, [2,1,3]));
        end


        function obj = run_parallel(obj, mtexHREBSD)
            obj.progressBar = waitbar(0, "Starting parallel pools...");
            switch obj.methodAnalysis
                case 'classic'
                    obj = obj.run_analysisParallel(mtexHREBSD);
                case 'hybrid'
                    obj = obj.run_analysisHybridParallel(mtexHREBSD);
                case 'dynamicSimulated'
                    delete(obj.progressBar);
                    obj = obj.run_analysisDynamicSimulatedParallel(mtexHREBSD);
                otherwise
                    disp("Analysis method " + string(obj.methodAnalysis) + " unknown.")
            end
        end


        obj = run_analysisDynamicSimulatedParallel(obj, mtexHREBSD)
        obj = run_analysisParallel(obj,mtexHREBSD)
        obj = run_grainSetRef(obj, mtexHREBSD, refPat)


        function obj = run_notParallel(obj, mtexHREBSD)
            obj.progressBar = waitbar(0, "Starting analysis...");
            switch obj.methodAnalysis
                case 'classic'
                    obj = obj.run_analysis(mtexHREBSD);
                otherwise
                    disp("Analysis method " + string(obj.methodAnalysis) + " unknown.")
            end
        end


        function obj = initialize(obj, mtexHREBSD)
            obj.version = mtexHREBSD.version;
            scanLength = mtexHREBSD.scan.scanLength;
            obj.beta = zeros([3,3,scanLength]);
            obj.F = zeros([3,3,scanLength]);
            obj.ft = mtexHREBSD.ft;
            obj.g = zeros([3,3,scanLength]);
            obj.fitMetrics = struct("SSE", zeros(1, scanLength), ...
                                    "R2", zeros(1, scanLength));
            obj.logMessage = {};
        end


        function update_progress(obj, iter, total, text)
            waitbar(iter/total, obj.progressBar, text);
        end


        function obj = run_grain(obj, mtexHREBSD)
            switch obj.methodAnalysis
                case 'setReference'
                    obj = obj.run_grainSetRef(mtexHREBSD);
                case 'classic'
                    obj = obj.run_grainClassic(mtexHREBSD);
                case 'hybrid'
                    obj = obj.run_grainHybrid(mtexHREBSD);
                case 'dynamicSimulated'
                    obj = obj.run_grainDynamic(mtexHREBSD);
            end
        end


        function obj = run_grainDynamic(obj, mtexHREBSD)
            tic
            disp('Dynamic simulated')
            numPoints = length(mtexHREBSD.ebsd);
%             numPoints = 10;
            obj.progressBar = waitbar(0, "Starting parallel pools...");
            mtexHREBSD_main.nUpdateWaitbar(numPoints, obj.progressBar)
            data = cell(numPoints, 1);
%             parpool(mtexHREBSD.analysis.numCores)
            D = parallel.pool.DataQueue;
            afterEach(D, @mtexHREBSD_main.nUpdateWaitbar);
            parfor j = 1:numPoints
                try
                    temp = mtexHREBSD;
                    testPat = temp.get_pattern(j);
                    calcF = dynamicSimulatedCalcF(testPat, temp);
                    data{j} = calcF;
                catch ME
                    disp(ME)
                end
                send(D, 1);
            end
            for i = 1:numPoints
                try
                    obj.beta(:,:,i) = data{i}.beta; 
                    obj.F(:,:,i) = data{i}.F; 
                    obj.g(:,:,i) = data{i}.g;
                    obj.fitMetrics.SSE(i) = data{i}.fit.metrics.SSE;
                    obj.fitMetrics.R2(i) = data{i}.fit.metrics.R2;
                catch ME
                    disp(ME)
                end
            end
            obj.completionTime = toc;
        end


        function obj = run_grainHybrid(obj, mtexHREBSD)
            disp('Hybrid')
            tic
            numPoints = length(mtexHREBSD.ebsd);
%             obj.globalIds = zeros(numPoints,1);
            localRefId = find(mtexHREBSD.ebsd.id == mtexHREBSD.refIds);
            refPat = mtexHREBSD.get_pattern(localRefId);
            obj.refF = dynamicSimulatedCalcF(refPat, mtexHREBSD,'ConvergencePlots',1);
            if ~isempty(obj.refF.logMessage)
                obj.logMessage{end+1} = obj.refF.logMessage;
            end
%             p = waitbar(0/numPoints, 'Starting hybrid analysis');
            for j = 1:numPoints
                try
                    testPat = mtexHREBSD.get_pattern(j);
%                     obj.globalIds(j) = testPat.scanIndex;
                    calcF = classicCalcF(refPat, testPat, mtexHREBSD);
%                     calcF.beta = calcF.beta - obj.refF.beta;
%                     calcF.F = calcF.beta + eye(3);
%                     calcF.beta = calcF.beta*inv(refF.beta);

                    calcF.F = calcF.F*obj.refF.F;
%                     calcF.F = obj.refF.F\calcF.F;
                    calcF.beta = calcF.F - eye(3);

%                     diffF = calcF.F + refF.F;
                    obj.beta(:,:,j) = calcF.beta;
                    obj.F(:,:,j) = calcF.F;
                    obj.g(:,:,j) = calcF.g;
                    obj.fitMetrics.SSE(j) = calcF.fit.metrics.SSE;
                    obj.fitMetrics.R2(j) = calcF.fit.metrics.R2;
                catch ME
                    obj.logMessage{end+1} = ME;
%                     disp(ME)
                end
%                 if ~mod(j, 10)
%                     time_j = toc;
%                     obj.update_waitbar(time_j, numPoints, j, p)
%                 end
            end
            time_j = toc;
%             obj.update_waitbar(time_j, numPoints, j, p)
        end
    end


    methods(Static)
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
                waitStr = "Completed points " + num2str(COUNT) + "/" + num2str(TOTAL) + ...
                    ", elapsed time " + num2str(h) + ":" + num2str(m) + ...
                    ":" + num2str(s);
                waitbar(p, H, waitStr)
            end
        end


        function update_waitbar(time, numPoints, j, p)
            [h,m,s] = hms(seconds(time));
            waitStr = "Points complete: " + num2str(j) + ...
                ", elapsed time: " + num2str(h) + ":" + num2str(m) + ":" + num2str(s);
            waitbar(j/numPoints, p, waitStr);
        end


        function methodCalcF = get_methodCalcF(mtexHREBSD)
            switch mtexHREBSD.analysis.calcFMethod
                case 'FDelta'
                    methodCalcF = 'calcFDelta';
                case 'classic'
                    methodCalcF = 'classicCalcF';
                otherwise
                    disp(['Calc F method unkown, defaulting to classic ' ...
                        'cross-correlation against real patterns']);
                    methodCalcF = 'classicCalcF';
            end
        end


        function methodAnalysis = get_methodAnalysis(mtexHREBSD)
            switch mtexHREBSD.analysis.analysisMethod
                case 'FDelta'
                    methodAnalysis = 'FDelta';
                case 'classic'
                    methodAnalysis = 'classic';
                case 'hybrid'
                    methodAnalysis = 'hybrid';
                case 'dynamicSimulated'
                    methodAnalysis = 'dynamicSimulated';
                otherwise
                    disp(['Analysis method unkown, defaulting to classic ' ...
                        'cross-correlation against real patterns']);
                    methodAnalysis = 'classic';
            end
        end
    end
end