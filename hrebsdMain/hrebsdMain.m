classdef hrebsdMain 
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        version
        beta
        F
        Fa
        Fc
        ft
        g
        fitMetrics
        logMessage
        args
        runTime
        peakHeights
    end


    methods
        function obj = hrebsdMain(mtexHREBSD, varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj.initialize(mtexHREBSD);
            obj.args = obj.get_args(varargin{:});
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


        function options = get_args(obj, varargin)
            p = obj.create_inputParser;
            parse(p, varargin{:});
            options = p.Results;
        end


        function obj = readAnalysis(obj, calcFChunk)
            for i = 1:size(calcFChunk,1)
                try
                    calcF = calcFChunk{l};
                    scanIndex = calcF.scanIndex;
                    obj.F(:,:,scanIndex) = calcF.F;
                    obj.g(:,:,scanIndex) = calcF.g;
                    obj.beta(:,:,scanIndex) = calcF.beta;
                    obj.fitMetrics.SSE(scanIndex) = calcF.fit.metrics.SSE;
                    obj.fitMetrics.R2(scanIndex) = calcF.fit.metrics.R2;
                catch me
                    obj.logMessage{end+1} = me;
                end
            end
        end
    end

    methods (Static)
        function p = create_inputParser
            p = inputParser;
            addParameter(p, 'Subset', 0);
        end

        function ebsd = ebsdSubset(mtexHREBSD, args)
            if args.Subset
                ebsd = mtexHREBSD.ebsd(args.Subset);
            else
                ebsd = mtexHREBSD.ebsd;
            end
        end

        function setup_parallel(mtexHEBSD)
            pool = gcp("nocreate");
            if isempty(pool)
                parpool(mtexHEBSD.analysis.numCores);
            end
        end
    end
end