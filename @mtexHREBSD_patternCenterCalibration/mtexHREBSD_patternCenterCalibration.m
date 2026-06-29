classdef mtexHREBSD_patternCenterCalibration < dynamicprops
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        gridIndicies
        method
        gridSpacing
        initialPC
        actualPC
        location
        elapsedTime
    end

    properties(Hidden)
        progressBar
        progressBarTitle
    end


    methods
        function obj = mtexHREBSD_patternCenterCalibration(mtexHREBSD, varargin)
            obj.progressBarTitle = inputname(1);
            obj = obj.initialize(varargin{:});
            if isempty(obj.gridIndicies)
                obj.gridSpacing = obj.get_gridSpacing(mtexHREBSD);
                obj.gridIndicies = obj.get_gridIndicies(mtexHREBSD);
            end
        end

%         [patternCenter, initialPatternCenter, location] = ...
%             PCCalParallel(obj, mtexHREBSD, kwargs)


        function obj = get_patternCenterCal(obj, mtexHREBSD, kwargs)
            arguments
                obj mtexHREBSD_patternCenterCalibration
                mtexHREBSD 
                kwargs.Verbose {mustBeMember(kwargs.Verbose, [0,1])} = 0
                kwargs.Simulators = []
            end

            disp("Starting PC Calibration...")
%             simulator = EBSPSim(mtexHREBSD, 'Ni');
            if ~isempty(kwargs.Simulators)
                calcf = mtexHREBSD.initializeCalcF( ...
                    'AnalysisType', 2, ...
                    'Simulators', kwargs.Simulators ...
                    );
            else
                calcf = mtexHREBSD.initializeCalcF( ...
                    'AnalysisType', 2 ...
                    );
            end
            x = mtexHREBSD.ebsd.prop.x(obj.gridIndicies);
            y = mtexHREBSD.ebsd.prop.y(obj.gridIndicies);
            obj.location = [x,y];
            if mtexHREBSD.analysis.numCores > 1
                outPCCal = mtexHREBSDMain(mtexHREBSD,...
                    'calcFStrain',calcf,...
                    'Subset',obj.gridIndicies, ...
                    'ParallelScheme',2, ...
                    'ChunkSize',1, ...
                    'Verbose', kwargs.Verbose, ...
                    'numCores', mtexHREBSD.analysis.numCores ...
                    );
            else
                outPCCal = mtexHREBSDMain(mtexHREBSD,...
                    'calcFStrain',calcf,...
                    'Subset',obj.gridIndicies ...
                    );
            end
            obj.initialPC = outPCCal.PCInitial(:,obj.gridIndicies)';
            obj.actualPC = outPCCal.PCCalibrated(:,obj.gridIndicies)';
                

%             if mtexHREBSD.analysis.numCores > 1
%                 [patternCenter, initialPatternCenter, loc] = ...
%                     obj.PCCalParallel(mtexHREBSD, kwargs);
%             else
%                 [patternCenter, initialPatternCenter, loc] = ...
%                     obj.get_patternCenterFDelta(mtexHREBSD);
%             end
            obj.elapsedTime = toc;
%             obj.initialPC = initialPatternCenter;
%             obj.actualPC = patternCenter;
%             obj.location = loc;
        end


        function patternCenterOffset = get_patternCenterOffset(obj, kwargs)
            arguments
                obj mtexHREBSD_patternCenterCalibration
                kwargs.FitType {mustBeMember(kwargs.FitType,...
                                        ["Normal", "Kernel"])} = "Normal"
            end
            patternCenterOffset = zeros(1,3);
            for i = 1:3
                offsetDist = rmoutliers(obj.actualPC(:,i) - obj.initialPC(:,i));
                if strcmp(kwargs.FitType, 'Normal')
                    pd = fitdist(offsetDist, 'Normal');
                    patternCenterOffset(i) = pd.mu;
                elseif strcmp(kwargs.FitType, 'Kernel')
                    pd = fitdist(offsetDist, 'Kernel');
                    patternCenterOffset(i) = pd.mean;
                end
            end
        end


        function [patternCenter, initialPatternCenter, location] = get_patternCenterFDeltaParallel(obj, mtexHREBSD)
            disp("Starting PC Calibration...")
            tic
            pool = gcp("nocreate");
            if isempty(pool)
                parpool(mtexHREBSD.analysis.numCores);
            end
            wb = waitbar(0, "Starting analysis", "Name",obj.progressBarTitle);
            patternCenter = zeros(length(obj.gridIndicies), 3);
            initialPatternCenter = zeros(length(obj.gridIndicies), 3);
            location = zeros(length(obj.gridIndicies), 2);
            testPats = cell(length(obj.gridIndicies), 1);
            for i = 1:length(obj.gridIndicies)
                waitbar(i/length(obj.gridIndicies), wb, "Getting patterns...")
                testPats{i} = mtexHREBSD.get_pattern(obj.gridIndicies(i));
                initialPatternCenter(i,:) = testPats{i}.patternCenter;
                location(i,1) = mtexHREBSD.ebsd.prop.x(obj.gridIndicies(i));
                location(i,2) = mtexHREBSD.ebsd.prop.y(obj.gridIndicies(i));
            end
            C = parallel.pool.Constant(mtexHREBSD);
            D = parallel.pool.DataQueue;
            waitbar(0, wb, "Starting pattern center calibration...")
            nUpdateWaitbar(length(obj.gridIndicies), wb)
            afterEach(D, @nUpdateWaitbar);
            tic  
            parfor i = 1:length(obj.gridIndicies)
                try   
                    calcf = FDeltaCalcF(testPats{i}, C.Value);
                catch me
                    disp(me)
                    calcf = struct('patternCenter', [0,0,0]);
                end
                patternCenter(i,:) = calcf.patternCenter;
                send(D,1)
            end
            toc
            disp("Finished PC Calibration!")
        end


        function [patternCenter, initialPatternCenter, location] = get_patternCenterFDelta(obj, mtexHREBSD)
            disp("Starting PC Calibration...")
            tic
            wb = waitbar(0, "Starting analysis", "Name",obj.progressBarTitle);
            patternCenter = zeros(length(obj.gridIndicies), 3);
            initialPatternCenter = zeros(length(obj.gridIndicies), 3);
            location = zeros(length(obj.gridIndicies), 2);
            D = parallel.pool.DataQueue;
            waitbar(0, wb, "Starting pattern center calibration...")
            nUpdateWaitbar(length(obj.gridIndicies), wb)
            afterEach(D, @nUpdateWaitbar);
            tic  
            for i = 1:length(obj.gridIndicies)
                testPat = mtexHREBSD.get_pattern(obj.gridIndicies(i));
                initialPatternCenter(i,:) = testPat.patternCenter;
                location(i,1) = mtexHREBSD.ebsd.prop.x(obj.gridIndicies(i));
                location(i,2) = mtexHREBSD.ebsd.prop.y(obj.gridIndicies(i));
                try
                    calcf = FDeltaCalcF(testPat, mtexHREBSD);
                catch me
                    disp(me)
                    calcf = struct('patternCenter', [0,0,0]);
                end
                patternCenter(i,:) = calcf.patternCenter;
                send(D,1)
            end
            toc
            disp("Finished PC Calibration!")
        end


        function gridIndicies = get_gridIndicies(obj, mtexHREBSD)
            xStop = round(mtexHREBSD.scan.Nx);
            yStop = round(mtexHREBSD.scan.Ny);
            xPoints = (0:obj.gridSpacing:xStop-1)*mtexHREBSD.scan.xStep;
            yPoints = (0:obj.gridSpacing:yStop-1)*mtexHREBSD.scan.yStep;
            gridIndicies = zeros(length(xPoints)*length(yPoints),1);
            i = 0;
            for x = xPoints
                for y = yPoints
                    ebsdX = mtexHREBSD.ebsd.prop.x;
                    ebsdY = mtexHREBSD.ebsd.prop.y;
                    ebsdX(abs(ebsdX) < 1E-5) = 0;
                    ebsdY(abs(ebsdY) < 1E-5) = 0;
                    i = i + 1;
                    xLogic = abs(abs(ebsdX) - x) < 1E-5;
                    yLogic = abs(abs(ebsdY) - y) < 1E-5;
%                     disp(i)
                    ind = find(xLogic & yLogic, 1);
                    gridIndicies(i) = ind;
                end
            end
        end


        function gridSpacing = get_gridSpacing(obj, mtexHREBSD)
            if obj.gridSpacing
                gridSpacing = obj.gridSpacing;
            else
                gridSpacing = round(mtexHREBSD.scan.Nx/10);
            end
        end


        
        function obj = initialize(obj, varargin)
            options = obj.parse_inputs(varargin{:});
            obj.method = options.method;
            obj.gridSpacing = options.gridSpacing;
            if isempty(options.indicies)
                obj.gridSpacing = options.gridSpacing;
            else
                obj.gridIndicies = options.indicies;
            end
        end


        function progressBar = start_waitbar(obj)
            waitString = "Starting pattern center calibration using " + ...
                          obj.method + "...";
            progressBar = waitbar(0, waitString);
        end


        function update_waitbar(obj, i)
            ind = obj.gridIndicies(i);
            waitString = "Pattern center calibration, current index: " + num2str(ind);
            waitbar(i/length(obj.gridIndicies), obj.progressBar, waitString);
        end


%         function pccalplots(obj, ebsd)
%             figure
%             plot(ebsd)
%             legend off
%             hold on
%             scatter( ...
%                 ebsd.prop.x(obj.gridIndicies), ...
%                 ebsd.prop.y(obj.gridIndicies), ...
%                 'k.');
%             for i = 1:length(obj.gridIndicies)
%                 x = ebsd.prop.x(obj.gridIndicies(i));
%                 y = ebsd.prop.y(obj.gridIndicies(i));
%                 text(x-1,y-0.5,num2str(i))
%             end
%         end


        function plot(obj, varargin)
            scatter3(obj.initialPC(:,1), obj.initialPC(:,3), obj.initialPC(:,2), '.r');
            hold on
            h = scatter3(obj.actualPC(:,1), obj.actualPC(:,3), obj.actualPC(:,2), 5, '.b');
            set(h, 'MarkerEdgeAlpha', 0.1, 'MarkerFaceAlpha', 0.1)
            pcRange = obj.get_patternCenterRange(obj.actualPC, 1);
            if nargin == 2
                shiftedPC = obj.initialPC + varargin{1};
                scatter3(shiftedPC(:,1), shiftedPC(:,3), shiftedPC(:,2), '.k');
                totalData = [obj.initialPC; pcRange; shiftedPC];
                axis equal
                xlabel('x*')
                ylabel('z*')
                zlabel('y*')
                legend("Uncalibrated", "F\Delta", "Calibrated")
            else
                totalData = [obj.initialPC; pcRange];
            end
            hold off
            axis equal
            xlim([min(totalData(:,1)), max(totalData(:,1))])
            ylim([min(totalData(:,3)), max(totalData(:,3))])
            zlim([min(totalData(:,2)), max(totalData(:,2))])
            if nargin == 2
                totalData2 = [pcRange; shiftedPC];
                figure
                scatter(obj.actualPC(:,1), obj.actualPC(:,2), '.b');
                hold on
                scatter(shiftedPC(:,1), shiftedPC(:,2), '.k');
                hold off
                axis equal
                grid on
                xlim([min(totalData2(:,1)), max(totalData2(:,1))])
                ylim([min(totalData2(:,2)), max(totalData2(:,2))])
            end
        end
    end


    methods(Static)
        function patternCenterRange = get_patternCenterRange(patternCenters, numStds)
            pcStd = std(patternCenters);
            pcMean = mean(patternCenters);
            logicArrays = cell(1, 3);
            for i = 1:3
                logicArrays{i} = patternCenters(:,i) > pcMean(i) - numStds*pcStd(i) & patternCenters(:,i) < pcMean(i) + numStds*pcStd(i);
            end
            indices = logicArrays{1} & logicArrays{2} & logicArrays{3};
            patternCenterRange = patternCenters(indices,:);
        end

        function options = parse_inputs(varargin)
            p = inputParser;
            addParameter(p,'method', 'F-Delta')
            addParameter(p,'gridSpacing', 0)
            addParameter(p,'indicies', [])
            parse(p, varargin{:})
            options = p.Results;
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
                [hr,m,s] = hms(seconds(t));
                waitStr = "Completed points " + num2str(COUNT) + ...
                    ", elapsed time " + num2str(hr) + ":" + num2str(m) + ...
                    ":" + num2str(s);
                waitbar(p, H, waitStr)
            end
        end
    end
end