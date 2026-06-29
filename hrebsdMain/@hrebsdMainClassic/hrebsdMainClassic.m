classdef hrebsdMainClassic < hrebsdMain
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here


    methods
        function obj = hrebsdMainClassic(mtexHREBSD, varargin)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@hrebsdMain(mtexHREBSD, varargin{:});
            obj = obj.runAnalysis(mtexHREBSD);
        end


        function obj = runAnalysis(obj, mtexHREBSD)
            timer = tic;
            if mtexHREBSD.analysis.numCores > 1
                obj.setup_parallel(mtexHREBSD)
                obj = obj.runParallel2(mtexHREBSD);
%                 obj = obj.runParallel3(mtexHREBSD);
            else 
                obj = obj.runStandard(mtexHREBSD);
            end
            obj.runTime = toc(timer);
        end

        
        obj = runParallel2(obj, mtexHREBSD);
        obj = runParallel3(obj, mtexHREBSD);
        obj = runStandard(obj, mtexHREBSD);
        obj = runStandardDislocation()
    end
end