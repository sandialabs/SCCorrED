classdef calcF1F2 < mtexHREBSD_calcF
    %calcF1F2 do the thing for the derivatives
    %   Detailed explanation goes here

    properties
        inds
    end

    methods
        function obj = calcF1F2(testPattern, mtexHREBSD, varargin)
            obj = obj@mtexHREBSD_calcF(mtexHREBSD);
            obj.inds = obj.get_xyIndicies(testPattern, mtexHREBSD);
            obj = obj.performAnalysis(testPattern, mtexHREBSD,varargin{:});
%             scanIndex = testPattern.scanIndex;
        end


        function obj = performAnalysis( ...
                obj, testPattern, mtexHREBSD, varargin ...
                )
            obj.beta = zeros([3,3,2]);
            obj.F = zeros([3,3,2]);
            obj.g = zeros([3,3,2]);
            obj.peakHeight = [0,0];
            obj.qs = zeros([3, mtexHREBSD.roi.numRois,2]);
            obj.rs = zeros([3, mtexHREBSD.roi.numRois,2]);
            for i = 1:length(obj.inds)
                pat = mtexHREBSD.get_pattern(obj.inds(i));
                calcF = classicCalcF( ...
                    testPattern, pat, mtexHREBSD, varargin{:} ...
                    );
                if i == 1
                    obj.A = zeros([size(calcF.A,1), size(calcF.A,2), 2]);
                    obj.b = zeros([size(calcF.b,1), size(calcF.b,2), 2]);
                end
                obj.A(:,:,i) = calcF.A;
                obj.b(:,:,i) = calcF.b;
                obj.F(:,:,i) = calcF.F;
                obj.beta(:,:,i) = calcF.beta;
                obj.g(:,:,i) = calcF.g;
                obj.peakHeight(i) = calcF.peakHeight;
                obj.qs(:,:,i) = calcF.qs;
                obj.rs(:,:,i) = calcF.rs;
            end
            obj.doRobustFit = calcF.doRobustFit;
        end
    end


    methods(Static)
        function inds = get_xyIndicies(pattern, mtexHREBSD)
            inds = [0,0];
            xstep = mtexHREBSD.scan.xStep;
            ystep = mtexHREBSD.scan.yStep;
            inds(1) = mtexHREBSD.findIndexXY( ...
                [pattern.xy(1) + xstep, pattern.xy(2)] ...
                );
            inds(2) = mtexHREBSD.findIndexXY( ...
                [pattern.xy(1), pattern.xy(2) + ystep] ...
                );
        end
    end
end