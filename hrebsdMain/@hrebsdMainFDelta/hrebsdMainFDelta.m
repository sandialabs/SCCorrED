classdef hrebsdMainFDelta < hrebsdMain
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        residual
        pcs
    end

    methods
        function obj = hrebsdMainFDelta(mtexHREBSD, varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            mtexHREBSD.analysis.assumptions = 'trace=0';
            obj = obj@hrebsdMain(mtexHREBSD, varargin{:});
            obj = obj.runAnalysis(mtexHREBSD);
        end

        function obj = runAnalysis(obj, mtexHREBSD)
            timer = tic;
            if mtexHREBSD.analysis.numCores > 1
                obj.setup_parallel(mtexHREBSD)
                obj = obj.runParallelFDelta(mtexHREBSD);
            end
            obj.runTime = toc(timer);
        end

        obj = runParallelFDelta(obj, mtexHREBSD)
    end
end