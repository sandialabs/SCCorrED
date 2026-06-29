classdef mtexHREBSD_dummy
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        analysis
        scan
        roi
        ft
        C
        firstPattern
        paths
    end

    methods
        function obj = mtexHREBSD_dummy(mtexHREBSD, dummyPattern, varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.scan = mtexHREBSD.scan;
            obj.analysis = mtexHREBSD.analysis;
            obj.firstPattern = dummyPattern;
            obj.roi = mtexHREBSD_roiSettings(obj.analysis, obj.firstPattern.image);
            obj.scan.pixelSize = obj.roi.pixelSize;
            obj.C = mtexHREBSD.C;
            obj.ft = mtexHREBSD.ft;
            obj.ft.m = obj.roi.pixelSize;
            obj.paths = mtexHREBSD.paths;
        end


        function fileOptions = get_fileOptions(obj, varargin)
            inds = obj.parse_inputs(varargin, "fileOptions");
            if ~isempty(inds)
                fileOptions = varargin{inds(2)};
            else 
                fileOptions = {};
            end
        end


        function analysisOptions = get_analysisOptions(obj, varargin)
            inds = obj.parse_inputs(varargin, "analysisOptions");
            if ~isempty(inds)
                analysisOptions = varargin{inds(2)};
            else 
                analysisOptions = {};
            end
        end
    end


    methods(Static)
        function out = parse_inputs(input, which_option)
            ind = find(strcmp(input, string(which_option)));
            out = [ind, ind+1];
        end
    end
end