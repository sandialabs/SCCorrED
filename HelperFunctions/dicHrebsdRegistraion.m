classdef dicHrebsdRegistraion < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numPoints   double
        imageRes    double % pixel/micron
        dicPoints   double
        ebsdPoints  double
        rotation    double
        offset      double
        type        
    end
    

    methods
        function obj = dicHrebsdRegistraion(numPoints, imageRes)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.numPoints = numPoints;
            obj.imageRes = imageRes;
        end


        function [X,Y] = transformPoints(obj,x,y)
            X = zeros(size(x));
            Y = zeros(size(y));
            for i = 1:length(X)
                tempCoord = [x(i), y(i), 0]';
                xform = obj.rotation*tempCoord + obj.offset;
                X(i) = xform(1);
                Y(i) = xform(2);
            end
        end


        function [X,Y] = transformDicCoords(obj, dicData)
            xdic = obj.dicPoints;
            xdic(:,3) = 1;
            xebsd = obj.ebsdPoints;
            xebsd(:,3) = 1;
            F = (xdic\xebsd);
            u = dicData.u.*obj.imageRes;
            v = dicData.v.*obj.imageRes;
            x = dicData.x.*obj.imageRes + u;
            y = dicData.y.*obj.imageRes + v;
            X = zeros(size(x));
            Y = zeros(size(y));
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                    tempCoord = [x(i,j), y(i,j), 1]';
                    tformCoord = F'*tempCoord;
                    X(i,j) = tformCoord(1);
                    Y(i,j) = tformCoord(2);
                end
            end
        end


%         function [X,Y] = transformDicCoords(obj, dicData)
%             u = dicData.u.*obj.imageRes;
%             v = dicData.v.*obj.imageRes;
%             x = dicData.x.*obj.imageRes + u;
%             y = dicData.y.*obj.imageRes + v;
%             X = zeros(size(x));
%             Y = zeros(size(y));
%             for i = 1:size(X,1)
%                 for j = 1:size(X,2)
%                     tempCoord = [x(i,j); y(i,j); 0];
%                     transformedCoord = obj.rotation*tempCoord + obj.offset;
%                     X(i,j) = transformedCoord(1);
%                     Y(i,j) = transformedCoord(2);
%                 end
%             end
%         end

        function obj = calcRegistratinCoefficents(obj)
            [obj.rotation, obj.offset] = fiducial_register( ...
                obj.dicPoints', obj.ebsdPoints');
        end

        
        function obj = findDicPoints(obj, image)
            f = figure;
            imshow(image)
            points = zeros(obj.numPoints, 2);
            for i = 1:obj.numPoints
                roi = drawpoint;
                points(i,:) = roi.Position*obj.imageRes;
            end
            obj.dicPoints = points;
            close(f)
        end


        function obj = findEbsdPoints(obj, ebsd)
            f = figure;
            plot(ebsd, ebsd.prop.iq, 'micronbar','off')
            colormap gray
            points = zeros(obj.numPoints, 2);
            for i = 1:obj.numPoints
                roi = drawpoint;
                points(i,:) = roi.Position;
            end
            obj.ebsdPoints = points;
            close(f)
        end


        function set.rotation(obj, val)
            obj.rotation = val;
        end


        function set.offset(obj, val)
            obj.offset = val;
        end


        function set.dicPoints(obj, rois)
            obj.dicPoints = rois;
        end


        function set.ebsdPoints(obj, rois)
            obj.ebsdPoints = rois;
        end


        function set.type(obj, type)
            obj.type = type;
        end


        function subsetIndicies = findDataSubset(obj, dic, hrebsd)
            [x,y] = obj.transformDicCoords(dic);
            area = obj.findArea(x,y);
            subsetIndicies = obj.findSubset(area, hrebsd);
        end


        function subsetIndicies = findDataSubsetFullGrains(obj, dic, hrebsd, minGrainSize)
            [x,y] = obj.transformDicCoords(dic);
            area = obj.findArea(x,y);
            subsetIndicies2 = obj.findSubset(area, hrebsd);
            ebsdTemp = hrebsd.ebsd(subsetIndicies2);
            ebsd = ebsdTemp(inpolygon(ebsdTemp.prop.x, ebsdTemp.prop.y, area(1,:), area(2,:)));
            grains = unique(ebsd.grainId);
            fullGrains = [];
            for i = 1:length(grains)
                ebsdGrain = hrebsd.ebsd(hrebsd.ebsd.grainId == grains(i));
                check = inpolygon(ebsdGrain.prop.x, ebsdGrain.prop.y, area(1,:), area(2,:));
                if all(check)
                    if hrebsd.grains.grainSize(grains(i)) > minGrainSize
                        fullGrains = [fullGrains, grains(i)];
                    end
                end
            end
            subsetIndicies = [];
            for j = 1:length(fullGrains)
                subsetIndicies = [subsetIndicies; hrebsd.ebsd.id(hrebsd.ebsd.grainId == fullGrains(j))];
            end
        end

    end
    methods(Static)
%         function subsetIndicies = findDataSubset(reg, dic, hrebsd)
%             [x,y] = reg.transformDicCoords(dic);
%             area = findArea(x,y);
%             subsetIndicies = findSubset(area, hrebsd);
%         end
        
        function area = findArea(x,y)
            xMin = min(x(:));
            xMax = max(x(:));
            yMin = min(y(:));
            yMax = max(y(:));
            xRange = xMax - xMin;
            yRange = yMax - yMin;
            xMin = xMin + 0.025*xRange;
            xMax = xMax - 0.025*xRange;
            yMin = yMin + 0.025*yRange;
            yMax = yMax - 0.025*yRange;
            area = [xMin, xMax, xMax, xMin;...
                    yMin, yMin, yMax, yMax];
        end
        
        function subsetIndicies = findSubset(area, hrebsd)
            ebsdArea = inpolygon(hrebsd.ebsd.prop.x, hrebsd.ebsd.prop.y, area(1,:), area(2,:));
            subsetIndicies = hrebsd.ebsd.id(ebsdArea);
        end
    end
end