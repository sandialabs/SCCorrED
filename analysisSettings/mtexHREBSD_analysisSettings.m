classdef mtexHREBSD_analysisSettings
    %mtexHREBSD_analysisSettings Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numCores = 1
        doStrain = 1
        doDislocationDensity = 0
        roiStyle = 'Annular'
        roiShape = 'Square'
        roiFilter = [2,50,1,1]
        roiSizePercent = 25
        numRois = 48
        HROIMMethod = 'Real'
        calcFMethod = 'classic'
        analysisMethod = 'classic'
        assumptions = 'free-surface'
        GNDMethod = 'Full'
        DDSMethod = 'Nye-Kroner'
        imageFilterType = 'standard'
        imagefilter = [9,90,0,0]
        iterationLimit = 10
        stepTolerance = 1.0E-4
        calcDerivatives = 1
        numRefIds = 10
        initialPatternCenter = 'naive'
        patternCenterCalibration = 'FDelta'
        guessF = eye(3,3)
        outputPath = pwd; 
        simulationType = 'Kinematic'
        Gradient = 0
        ForceNaivePC = 0;
    end

    methods
        function obj = mtexHREBSD_analysisSettings(varargin)
            if nargin > 0
                input = varargin{:};
            else
                input = {};
            end
            props = properties(obj);
            for i = 1:length(props)
                obj = obj.update_var(string(props(i)), input);
            end
        end

        function obj = update_var(obj, prop, options)
            ind = find(strcmp(string(options), prop));
            if ~isempty(ind)
                obj.(prop) = options{ind+1};
            end
        end
    end
end