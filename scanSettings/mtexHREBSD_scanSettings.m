classdef (Abstract) mtexHREBSD_scanSettings
    %mtexHREBSD_scanSettings Summary of this class goes here
    %   Detailed explanation goes here

    properties
        cameraElevation = 10
        scanfile
        imagefile
        patternCenter = [0.5, 0.5, 0.5]
        xStep
        yStep
        Nx = 0
        Ny = 0
        scanLength
        workingDistance
        sampleTilt = 70
        KV = 20
        delta = 25
        pixelSize
        phosphorSize = 30000
        material = 'nickel'
    end

    methods
        function obj = mtexHREBSD_scanSettings(options)
            %mtexHREBSD_scanSettings Construct an instance of this class
            %   Detailed explanation goes here
            obj.scanfile = options.scanfile;
            obj.imagefile = options.imagefile;
            [~,~,ext] = fileparts(obj.scanfile);
            if any(strcmp(ext, {'.ctf', '.h5oina'}))
                obj.phosphorSize = 28800;
            else
                obj.phosphorSize = 30000;
            end
        end


        function obj = update_delta(obj, roiSettings)
            [~,~,ext] = fileparts(obj.scanfile);
            if any(strcmp(ext, {'.ctf', '.h5oina'}))
                obj.delta = 28800/(roiSettings.pixelSize/1024*1244);
            else
                obj.delta = 30000/roiSettings.pixelSize;
            end
        end


        function obj = update_pixelSize(obj, roiSettings)
            obj.pixelSize = roiSettings.pixelSize;
        end
    end
end