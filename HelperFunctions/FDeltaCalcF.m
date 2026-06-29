classdef FDeltaCalcF < mtexHREBSD_calcF 
    % FDeltaCalcF Summary of this class goes here
    %   Variable inputs
    %       Verbose displays, if 1 details for each iteration. Default 0
    %       ConvergencePlots, if 1 plots the convergence of iterations. 
    %           Default 0
    %       TroubleshootingPlots, if 1 will plot each iteration's reference
    %       pattern with shifts superimposed. Default 0
    %       Tolerance, tolerance for convergence. Default 1E-5

    properties
        patternCenter
        refRotation
    end


    properties(Hidden)
        options
        patternCenterAllIters
    end


    methods
        function obj = FDeltaCalcF(testPat, mtexHREBSD, varargin)
            %UNTITLED9 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@mtexHREBSD_calcF(mtexHREBSD);
            obj.options = obj.parseOptions(varargin{:});
            [refPat, obj.F, obj.patternCenterAllIters] = FDelta( ...
                testPat, mtexHREBSD, obj.options ...
                );
            obj.patternCenter = refPat.patternCenter;
%             obj.F = newF;
            obj.beta = obj.F - eye(3);
            obj.g = refPat.g;
            obj.refRotation = refPat.rotations;
            obj.C = testPat.C;
            tempPat = obj.makeTempPat(refPat, testPat);
            calcFTemp = classicCalcF(refPat, tempPat, mtexHREBSD);
            obj.peakHeight = calcFTemp.peakHeight;
            obj.rs = calcFTemp.rs;
            obj.qs = calcFTemp.qs;
            obj.fit.metrics.SSE = calcFTemp.fit.metrics.SSE;
            obj.fit.metrics.R2 = calcFTemp.fit.metrics.R2;
        end
    end

    methods(Static)
        function tempPat = makeTempPat(refPat, testPat)
            tempPat = testPat;
            tempPat.patternCenter = refPat.patternCenter;
            tempPat.rotations = refPat.rotations;
            tempPat.g = refPat.g;
        end


        function options = parseOptions(varargin)
            p = inputParser;
            addParameter(p, 'ConvergenceData', 0)
            addParameter(p, 'Verbose', 0)
            addParameter(p, 'ConvergencePlots', 0)
            addParameter(p, 'TroubleshootingPlots', 0)
            addParameter(p, 'Tolerance', 1E-5)
            addParameter(p, 'doRobustFit', 0)
            parse(p, varargin{:});
            options = p.Results;
        end
    end
end