classdef mtexHREBSD_scanSettings_h5oina < mtexHREBSD_scanSettings
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        detectorOrientation
    end

    methods
        function obj = mtexHREBSD_scanSettings_h5oina(options)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@mtexHREBSD_scanSettings(options);
            fname = options.scanfile;
            h = '/1/EBSD/Header/';
            obj.xStep = h5read(fname, [h, 'X Step']);
            obj.yStep = h5read(fname, [h, 'Y Step']);
            obj.Nx = h5read(fname, [h, 'X Cells']);
            obj.Ny = h5read(fname, [h, 'Y Cells']);
            obj.scanLength = obj.Nx * obj.Ny;
            obj.workingDistance = h5read( ...
                fname, [h, 'Working Distance'] ...
                );
            obj.sampleTilt = h5read(fname, [h, 'Tilt Angle']) * 180/pi;
            obj.KV = h5read(fname, [h, 'Beam Voltage']);
            obj.detectorOrientation = h5read( ...
                fname, [h, 'Detector Orientation Euler'] ...
                );
            obj.cameraElevation = obj.detectorOrientation(2)*180/pi - 90;
        end
    end
end