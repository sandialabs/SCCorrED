classdef mtexHREBS_fit
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        metrics struct;
        p
    end

    methods
        function obj = mtexHREBS_fit(varargin)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.p = obj.create_inputParser;
        end

%         function metrics = initialize(obj)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end

    methods (Static)
        function p = create_inputParser
            p = inputParser;
            addParameter(p, 'SSE', 0);
            addParameter(p, 'R2', 0);
        end
    end
end