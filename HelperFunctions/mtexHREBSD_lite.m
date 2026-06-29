classdef mtexHREBSD_lite
    %mtexHREBSD_lite used as an input to parallel processing methods to
    %reduce I/O of sending complete ebsd/grain data

    properties
        analysis
        scan
        roi
        ft
        C
        patternCenterOffset
        patterns
    end

    methods
        function obj = mtexHREBSD_lite(mtexHREBSD)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj.analysis = mtexHREBSD.analysis;
            obj.scan = mtexHREBSD.scan;
            obj.roi = mtexHREBSD.roi;
            obj.ft = mtexHREBSD.ft;
            obj.C = mtexHREBSD.C;
            obj.patternCenterOffset = mtexHREBSD.patternCenterOffset;
            obj.patterns = mtexHREBSD.patterns;
        end
    end
end