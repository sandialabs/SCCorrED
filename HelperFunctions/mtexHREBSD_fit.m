classdef mtexHREBSD_fit
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        metrics struct;
    end

    methods
        function obj = mtexHREBSD_fit(varargin)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.metrics = obj.initialize(varargin{:});
        end


        function obj = calc_metrics(obj, b, bFit)
            metricNames = fieldnames(obj.metrics);
            for i = 1:length(metricNames)
                fname = "mtexHREBSD_fit.get_" + metricNames{i};
                obj.metrics.(metricNames{i}) = feval(fname, b, bFit);
            end
        end
    end

    methods(Static)
        function R2 = get_R2(b, bFit)
            residual = ((b - bFit).^2);
            ssr = sum(residual);
            squares = (b - mean(b)).^2;
            tsm = sum(squares);
            R2 = 1 - ssr/tsm;
        end


        function SSE = get_SSE(b, bFit)
            residual = ((b - bFit).^2);
            SSE = sum(residual).^0.5;
        end


        function metrics = initialize(varargin)
            metrics = struct;
            for i = 1:nargin
                metrics.(varargin{i}) = [];
            end
        end
    end
end