classdef mtexHREBSD_settings
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        scan
    end

    methods
        function obj = mtexHREBSD_settings(options)
           obj.scan = mtexHREBSD_scanSettings(options);
        end
    end
end