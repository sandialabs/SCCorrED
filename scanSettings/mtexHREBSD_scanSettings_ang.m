classdef mtexHREBSD_scanSettings_ang < mtexHREBSD_scanSettings
    %mtexHREBSD_scanSettings_ang Summary of this class goes here
    %   Detailed explanation goes here

    properties
        grid
    end

    methods
        function obj = mtexHREBSD_scanSettings_ang(options)
            %mtexHREBSD_scanSettings_ang Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@mtexHREBSD_scanSettings(options);
            obj = obj.get_propertiesFromFile(12);
            obj.scanLength = obj.Nx * obj.Ny;
        end


        function obj = get_propertiesFromFile(obj, numItems)
            expression = '[xyz]-star|SampleTiltAngle|CameraElevationAngle|GRID|XSTEP|YSTEP|NCOLS_ODD|NCOLS_EVEN|NROWS|WorkingDistance';
            fid = fopen(obj.scanfile);
            propsFound = 0;
            for i = 1:100
                line = fgetl(fid);
                match_str = string(regexp(line, expression, 'match'));
                switch match_str
                    case 'x-star'
                        obj.patternCenter(1) = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'y-star'
                        obj.patternCenter(2) = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'z-star'
                        obj.patternCenter(3) = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'SampleTiltAngle'
                        obj.sampleTilt = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'CameraElevationAngle'
                        obj.cameraElevation = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'GRID'
                        obj.grid = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'NCOLS_ODD'
                        ncols_odd = obj.add_property(line);
                        propsFound = propsFound + 1;
                        if ncols_odd > obj.Nx
                            obj.Nx = ncols_odd;
                        end  
                    case 'NCOLS_EVEN'
                        ncols_even = obj.add_property(line);
                        propsFound = propsFound + 1;
                        if ncols_even > obj.Nx
                            obj.Nx = ncols_even;
                        end 
                    case 'NROWS'
                        obj.Ny = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'XSTEP'
                        obj.xStep = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'YSTEP'
                        obj.yStep = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'WorkingDistance'
                        obj.workingDistance = obj.add_property(line);
                        propsFound = propsFound + 1;
                end
                if propsFound >= numItems
                    break
                end
            end
        end
    end


    methods (Static)
        function prop = add_property(line)
            match_locations = regexp(line, ' ');
            ind = match_locations(end) + 1;
            prop = str2double(string(line(ind:end)));
            if isnan(prop)
                prop = string(line(ind:end));
            end
        end
    end
end