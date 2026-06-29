classdef mtexHREBSD_roiSettings
    %mtexHREBSD_roiSettings Summary of this class goes here
    %   Detailed explanation goes here

    properties
        roiStyle
        roiShape
        roiFilter
        roiSizePercent
        roiSize
        numRois
        pixelSize
        centers
        custfilt
        windowfunc
    end

    methods
        function obj = mtexHREBSD_roiSettings(analysisSettings, firstImage)
            %mtexHREBSD_roiSettings an instance of this class
            %   Detailed explanation goes here
            obj.roiStyle = analysisSettings.roiStyle;
            obj.roiShape = analysisSettings.roiShape;
            obj.roiFilter = analysisSettings.roiFilter;
            obj.roiSizePercent = analysisSettings.roiSizePercent;
            obj.numRois = analysisSettings.numRois;
            obj.pixelSize = min(size(firstImage));
            obj.roiSize = round(obj.pixelSize * analysisSettings.roiSizePercent/100);
            [obj.custfilt, obj.windowfunc] = obj.build_filters(analysisSettings);
            obj.centers = obj.get_centers(analysisSettings);
        end   

        
        function centers = get_centers(obj, analysisSettings)
            switch analysisSettings.roiStyle
                case 'Annular'
                    centers = obj.get_centersAnnular(analysisSettings);
                otherwise
                    warning('mtexHREBSD:invalid_roiStyle', ['roiStyle must ' ...
                        'be ''Annular'' or *list other roiStyles. ' ...
                        'Defaulting to Annular roiStyle.']);
                    centers = obj.get_centersAnnular(analysisSettings);
            end
        end
        

        function [custfilt, windowfunc] = build_filters(obj, analysisSettings)
            lowerrad = analysisSettings.roiFilter(1);
            upperrad = analysisSettings.roiFilter(2);
            L = obj.roiSize + 1;
            xc = round(L/2);
            yc = round(L/2);
            filt = zeros(obj.roiSize, obj.roiSize);
            i = 1:obj.roiSize;
            j = 1:obj.roiSize;
            IJ = meshgrid(i,j);
            dist = sqrt((IJ-ones(size(IJ)).*xc).^2+(IJ'-ones(size(IJ)).*yc).^2);
            filt(dist<lowerrad | dist>upperrad) = 1;
            if analysisSettings.roiFilter(4) == 1
                filt(dist>upperrad & dist<upperrad+13)=erf((dist(dist>upperrad & dist<upperrad+13)-upperrad)/13*pi);
            end
            if analysisSettings.roiFilter(3) == 1
                filt(dist<lowerrad & dist>lowerrad-13)=erf(-(dist(dist<lowerrad & dist>lowerrad-13)-lowerrad)/13*pi);
            end
            filt = 1-filt;
            custfilt = fftshift((filt));
            xc = L/2;
            yc = L/2;
            windowfunc = cos((IJ-ones(size(IJ)).*xc)*pi/obj.roiSize)...
                .* cos((IJ'-ones(size(IJ)).*yc)*pi/obj.roiSize);
        end


        function centers = get_centersAnnular(obj, analysisSettings)
            centers = zeros(2, analysisSettings.numRois);
            centers(:,1) = obj.pixelSize/2;
            angSpacing = 2*pi / (analysisSettings.numRois - 1);
            radius = floor((obj.pixelSize-0-obj.roiSize)/3);
            i = 1:(analysisSettings.numRois-1);
            dx = radius*cos((i-1)*angSpacing);
            dy = radius*sin((i-1)*angSpacing);
            centers(1, 2:end) = obj.pixelSize/2 + dx;
            centers(2, 2:end) = obj.pixelSize/2 + dy;
%             A = repmat([0.5;0.5], 1, obj.numRois);
%             centers = round(centers + A) - A;
        end


        function [roi, rrange, crange] = get_squareRoi(obj, pattern, roiInd)
            center2 = (obj.centers(:,roiInd));
            center = [center2(2), center2(1)];
            h = (obj.roiSize/2);
            rrange = round(center(1)-h):round(center(1)+h)-1;
            crange = round(center(2)-h):round(center(2)+h)-1;
            roi = pattern.image(rrange, crange);
%             roi = roi - mean(roi);
        end


        function [roi, rrange, crange] = get_squareRoiShifted(obj, pattern, roiInd, shift)
            center2 = (obj.centers(:,roiInd) + shift);
            center = [center2(2), center2(1)];
            h = (obj.roiSize/2);
            rrange = round(center(1)-h):round(center(1)+h)-1;
            crange = round(center(2)-h):round(center(2)+h)-1;
%             obj.check_range(rrange)
%             obj.check_range(crange)
            roi = pattern.image(rrange, crange);
%             roi = roi - mean(roi);
        end


        function check_range(obj, range)
            minRange = min(range);
            maxRange = max(range);
            if minRange < 0 
                disp("min too low!")
            elseif maxRange > obj.pixelSize
                disp("max too high!")
            end
        end
        

        function plot_centers(obj)
            figure
            scatter(obj.centers(1,:), obj.centers(2,:))
            axis equal
        end
    end
end