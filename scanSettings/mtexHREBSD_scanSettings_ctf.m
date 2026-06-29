classdef mtexHREBSD_scanSettings_ctf < mtexHREBSD_scanSettings
    %mtexHREBSD_scanSettings_ctf Summary of this class goes here
    %   Detailed explanation goes here
    properties
        detectorOrientation
    end


    methods
        function obj = mtexHREBSD_scanSettings_ctf(options)
            %mtexHREBSD_scanSettings_ctf Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@mtexHREBSD_scanSettings(options);
            obj = obj.get_propertiesFromFile(5);
            obj.scanLength = obj.Nx * obj.Ny;
        end
        

        function obj = get_propertiesFromFile(obj, numItems)
            expression = 'XCells|YCells|XStep|YStep|KV';
            fid = fopen(obj.scanfile);
            propsFound = 0;
            for i = 1:100
                line = fgetl(fid);
                match_str = string(regexp(line, expression, 'match'));
                switch match_str
                    case 'XCells'
                        obj.Nx = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'YCells'
                        obj.Ny = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'XStep'
                        obj.xStep = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'YStep'
                        obj.yStep = obj.add_property(line);
                        propsFound = propsFound + 1;
                    case 'KV'
                        data_struct = obj.split_props_line(line);
                        obj.KV = data_struct.KV;
                        obj.sampleTilt = data_struct.TiltAngle;
                        obj.workingDistance = data_struct.WorkingDistance;
                        obj.detectorOrientation = [...
                            data_struct.DetectorOrientationE1;...
                            data_struct.DetectorOrientationE2;...
                            data_struct.DetectorOrientationE3]*pi/180;
                        propsFound = propsFound + 1;
                end
                if propsFound >= numItems
                    break
                end
            end
        end
        
    end


    methods(Static)    
        function data_struct = split_props_line(line)
            line_string = string(line);
            split_string = split(line_string, '!');
            data_string = split_string(2);
            data_array = strsplit(data_string);
            data = reshape(data_array(2:end), [2, length(data_array(2:end))/2]);
            data_struct = struct();
            for i = 1:length(data)
                data_struct.(data(1,i)) = str2double(data(2,i));
            end
        end


        function prop = add_property(line)
            match_locations = regexp(line, '\t');
            if length(match_locations) > 1
                ind = match_locations(end) + 1;
            else
                ind = match_locations + 1;
            end
            prop = str2double(string(line(ind:end)));
            if isnan(prop)
                prop = string(line(ind:end));
            end
        end
    end
end