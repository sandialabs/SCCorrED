classdef hrebsdMainSimulated < hrebsdMain
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function obj = hrebsdMainSimulated(mtexHREBSD, varargin)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@hrebsdMain(mtexHREBSD, varargin{:});
            obj = obj.runAnalysis(mtexHREBSD);
        end

        function obj = runAnalysis(obj, mtexHREBSD)
            timer = tic;
            if mtexHREBSD.analysis.numCores > 1
                obj.setup_parallel(mtexHREBSD)
                obj = obj.runParallel(mtexHREBSD);
            else 
                obj = obj.runStandard(mtexHREBSD);
            end
            obj.runTime = toc(timer);
        end

        obj = runParallel(obj, mtexHREBSD);
        obj = runStandard(obj, mtexHREBSD);
    end
end