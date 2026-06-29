classdef dynamicSimulatedCalcF < mtexHREBSD_calcF
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        scanIndex
        options
        betaVec
        converged
        finalCalcF
        intialCalcF
    end

    methods
        function obj = dynamicSimulatedCalcF(testPattern, mtexHREBSD, varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here 
            obj = obj@mtexHREBSD_calcF(mtexHREBSD);
%             obj.g = testPattern.g;
            obj.options = obj.get_options(varargin{:});
            obj.doRobustFit = obj.options.doRobustFit;
            obj.F_guess = eye(3);
            obj.scanIndex = testPattern.scanIndex;
%             refPat = testPattern;
%             refPat.rotations(1) = refPat.rotations(1)+pi/40;
%             refPat.g = refPat.get_g_matrix(refPat.rotations(1), refPat.rotations(2), refPat.rotations(3));
%             refPat = refPat.get_simulatedPattern(mtexHREBSD);
%             obj = obj.performAnalysis(refPat, testPattern, mtexHREBSD);

            refPattern = testPattern.get_simulatedPattern;
            obj = obj.performAnalysis(refPattern, testPattern, mtexHREBSD);
        end


        function obj = performAnalysis(obj, refPat, testPat, mtexHREBSD)
            calcF = classicCalcF(refPat, testPat, mtexHREBSD);
            obj.intialCalcF = calcF;
            obj.F = calcF.F;
            normF0 = norm(calcF.F);
            if obj.options.Verbose
                disp("Iteration 0:")
                disp("norm(F) = " +num2str(normF0))
                disp("F =")
                disp(calcF.F)
                disp("g =")
                disp(calcF.g)
                disp(" ")
            end
            if obj.options.ConvergencePlots
                obj.betaVec = zeros(3,3,mtexHREBSD.analysis.iterationLimit);
            end
            if obj.options.TroubleshootingPlots
                ts = troubleshootingFunctions;
            end
            for i = 1:mtexHREBSD.analysis.iterationLimit
                refPat = obj.update_refPat(refPat, calcF, mtexHREBSD);
                testPat = obj.update_testPat(refPat, testPat);
                % Are we sure we don't want to update the test pattern
                % rotations???
                calcF = classicCalcF(refPat, testPat, mtexHREBSD);
                normF1 = norm(calcF.F);
                if obj.options.Verbose
                    disp("Iteration "+num2str(i)+":")
                    disp("norm(F) = " +num2str(normF1))
                    disp("F =")
                    disp(calcF.F)
                    disp("g =")
                    disp(calcF.g)
%                     disp(refPat.g)
                    disp(" ")
                end
                if obj.options.TroubleshootingPlots
                    qs = mtexHREBSD.ft.phosphorFrame2imageVec(calcF.qs);
%                     ts.plot_shift(refPat, mtexHREBSD, 5*qs);
%                     qs = mtexHREBSD.ft.phosphorFrame2imageVec(calcF.qs);
%             ts.plot_shift(refPat, mtexHREBSD, 5*qs);
                    qFit = mtexHREBSD.ft.phosphorFrame2imageVec(calcF.qFit);            
                    ts.plot_twoShifts(refPat, mtexHREBSD, qs, qFit,2)
                end
                obj.converged = obj.check_convergence(normF0, normF1, obj.options.Tolerance);
%                 converged = obj.check_convergence2(calcF, i, obj.options.Tolerance);
                if obj.converged
%                     obj.converged = 1;
                    obj.F = calcF.F;
                    obj.g = refPat.g;
                    obj.beta = calcF.beta;
                    obj.fit = calcF.fit;
                    obj.betaVec = obj.update_betaVec(i, calcF, obj.converged);
                    obj.peakHeight = calcF.peakHeight;
                    break
                else
                    normF0 = normF1;
                    obj.betaVec = obj.update_betaVec(i, calcF);
                end
            end
            if ~obj.converged
%                 obj.converged = 0;
                obj.F = calcF.F;
%                 obj.F_sample = calcF.F_sample;
%                 obj.F_crystal = calcF.F_crystal;
                obj.beta = calcF.beta;
                obj.fit = calcF.fit;
                obj.g = refPat.g;
                obj.qs = calcF.qs;
                obj.rs = calcF.rs;
                obj.logMessage = "Failed to converge using dynamic " + ...
                    "simulation methods for scan index " + ...
                    num2str(refPat.scanIndex);
            end
            if obj.options.ConvergencePlots
                obj.plot_betaVecConvergence;
            end
            obj.finalCalcF = calcF;
        end


        function options = get_options(obj, varargin)
            p = obj.create_inputParser;
            parse(p, varargin{:});
            options = p.Results;
        end


        function betaVec = update_betaVec(obj, i, calcF, varargin)
            betaVec = obj.betaVec;
            if nargin == 3
                betaVec(:,:,i) = calcF.beta;
            else
                betaVec(:,:,i) = calcF.beta;
                betaVec = betaVec(:,:,1:i);
            end
        end


        function plot_betaVecConvergence(obj)
            figure
            numIters = size(obj.betaVec, 3);
            labels = cell(9,1);
            hold on
            for i = 1:3
                for j = 1:3
                    toPlot = reshape(obj.betaVec(i,j,:), [numIters,1]);
                    plot(1:numIters, toPlot);
                    labels{(i-1)*3+j} = "\beta_{"+num2str(i) + num2str(j)+"}";
                end
            end
            grid on
            legend(labels, 'Location', 'eastoutside');
        end
        

        function check = check_convergence2(obj, calcF, i, tolerance)
            diffBeta = abs(obj.betaVec(i) - calcF.beta);
            diffBetaVIP = [diffBeta(1,1), diffBeta(1,2), diffBeta(1,3),...
                diffBeta(2,1), diffBeta(2,2), diffBeta(3,3)];
            if all(diffBetaVIP < tolerance)
                check = 1;
            else
                check = 0;
            end
        end
    end

    methods(Static)
        function plot_FVec(FVec)
            iterations = size(FVec, 3);
            figure
            hold on
            for i=1:3
                for j = 1:3
                    toPlot = squeeze(FVec(i,j,:));
                    disp(size(toPlot(:,1)))
                    plot(iterations, toPlot(:,1));
                end
            end
        end


        function check = check_convergence(normF0, normF1, tolerance)
            diff = abs(normF1 - normF0)/normF1;
            if diff < tolerance
                check = 1;
            else
                check = 0;
            end
        end


        function newTestPat = update_testPat(refPat, testPat)
            newTestPat = testPat;
            newTestPat.rotations = refPat.rotations;
            newTestPat.g = refPat.g;
        end


        function newRefPat = update_refPat(refPat, calcF, mtexHREBSD)
            copyPat = refPat;
            Qps = mtexHREBSD.ft.Qp2s;
            g = refPat.g;
            [rr, ~] = poldec(g*Qps*calcF.F*Qps'*g');
            copyPat.g = rr'*refPat.g;
            [copyPat.rotations(1),copyPat.rotations(2),copyPat.rotations(3)] = gmat2euler(copyPat.g);
            newRefPat = copyPat.get_simulatedPattern;
        end


        function p = create_inputParser
            p = inputParser;
            addParameter(p, 'Verbose', 0)
            addParameter(p, 'ConvergencePlots', 0)
            addParameter(p, 'TroubleshootingPlots', 0)
            addParameter(p, 'Tolerance', 1E-5)
            addParameter(p, 'doRobustFit', 0)
%             addParameter()
        end
    end
end