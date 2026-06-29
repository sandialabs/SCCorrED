classdef viableRefIds
    %viableRefIds Summary of this class goes here
    %   Detailed explanation goes here

    properties
        globalIds
        localIds
        firstLocal
        firstGlobal
    end
    properties(Hidden = true)
        numIds = 10;
        centroid
    end

    methods
        function obj = viableRefIds(grainEBSD, varargin)
            %viableRefIds Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 1
                obj.numIds = varargin{1};
            end
            obj = obj.check_numIds(grainEBSD);
            obj.localIds = obj.get_localIds(grainEBSD);
            obj.centroid = obj.get_centroid(grainEBSD);
            obj.firstLocal = obj.get_refinedFirstPassLocal(grainEBSD);
            obj.firstGlobal = obj.convert_local2global(grainEBSD, obj.firstLocal);
            obj.globalIds = obj.convert_local2global(grainEBSD, obj.localIds);
        end
        

        function obj = check_numIds(obj, ebsd)
            if length(ebsd) < obj.numIds
                obj.numIds = length(ebsd);
            end
        end


        function firstLocal = get_refinedFirstPassLocal(obj, grainEBSD)
            dist = obj.get_distFromCentroid(grainEBSD, obj.localIds);
            A = [dist, obj.localIds];
            B = sortrows(A, 1);
            firstLocal = B(1,2);
        end


        function localIds = get_localIds(obj, grainEBSD)
            localIds = zeros(obj.numIds,1);
            iq = grainEBSD.prop.iq;
            ids = 1:length(iq);
            A = [iq, ids'];
            B = sortrows(A, 1, "descend");
            localIds(1:obj.numIds) = B(1:obj.numIds, 2);
        end


        function dist = get_distFromCentroid(obj, grainEBSD, localFirstPass)
            diffx = obj.centroid(1) - grainEBSD.prop.x(localFirstPass);
            diffy = obj.centroid(2) - grainEBSD.prop.y(localFirstPass);
            dist = sqrt(diffx.^2 + diffy.^2);
        end
    end


    methods(Static)
        function centroid = get_centroid(grainEBSD)
            xc = mean(grainEBSD.prop.x);
            yc = mean(grainEBSD.prop.y);
            centroid = [xc, yc];
        end


        function globalId = convert_local2global(grainEBSD, localId)
            globalId = grainEBSD.id(localId);
        end
    end
end