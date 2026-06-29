classdef troubleshootingFunctions
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here

    methods(Static)
        function plot_PCCalPoints(pccal, mtexHREBSD)
            x = mtexHREBSD.ebsd.prop.x(pccal.grainIndicies);
            y = mtexHREBSD.ebsd.prop.y(pccal.grainIndicies);
            plot(mtexHREBSD.grains)
            hold on
            scatter(x,y,'r.')
            hold off
        end


        function plot_rois(pattern, mtexHREBSD)
            h = mtexHREBSD.roi.roiSize/2;
            canvas = zeros(size(pattern.image));
            for i = 1:mtexHREBSD.roi.numRois
                roi = mtexHREBSD.roi.get_squareRoi(pattern, i);
                center = mtexHREBSD.roi.centers(:,i);
                rangex = round(center(1)-h:center(1)+h-1);
                rangey = round(center(2)-h:center(2)+h-1);
                canvas(rangex, rangey) = roi;
                imshow(uint8(canvas))
            end
        end


        function plot_shift(pattern, mtexHREBSD, qs)
%             figure
            centers = mtexHREBSD.roi.centers;
            pattern.show_pattern;
            hold on
%             scatter(centers(2,:), centers(1,:),50, 'b.');
%             scatter(centers(2,:)+qs(2,:), centers(1,:)+qs(1,:),50, 'r.');
%             for i = 1:length(centers)
            q = quiver(centers(2,:), centers(1,:), qs(1,:), qs(2,:),'r', 'LineWidth', 1);
            q.AutoScale = 'off';
            hold off
        end


        function plot_twoShifts(pattern, mtexHREBSD, q1, q2, mult)
            if ~exist('mult', 'var')
                mult = 1;
            end
            centers = mtexHREBSD.roi.centers;
            pattern.show_pattern
%             figure
            hold on
            q1 = quiver(centers(1,:), centers(2,:), mult*q1(1,:), mult*q1(2,:), 'g');
            q1.AutoScale = 'off';
            q2 = quiver(centers(1,:), centers(2,:), mult*q2(1,:), mult*q2(2,:), 'r');
            q2.AutoScale = 'off';
            hold off
            axis equal
        end
            


        function plot_qAndFit(pattern, mtexHREBSD, qs, qFit)
            centers = mtexHREBSD.roi.centers;
            pattern.show_pattern;
            figure
            hold on
%             scatter(centers(2,:), centers(1,:),50, 'b.');
%             scatter(centers(2,:)+qs(2,:), centers(1,:)+qs(1,:),50, 'r.');
%             for i = 1:length(centers)
            q = quiver(centers(2,:), centers(1,:), qs(2,:), qs(1,:), 'LineWidth', 1, 'Color','b');
            q2 = quiver(centers(2,:), centers(1,:), qFit(2,:), qFit(1,:), 'LineWidth', 1,'Color','r');
            legend('q', 'q_{fit}')
%             q.AutoScale = 'off';
%             q2.AutoScale = 'off';
            hold off
        end


        function plot_shift2(refpattern, pattern, mtexHREBSD, qs)
%             figure
            centers = mtexHREBSD.roi.centers;
            I1 = uint8(refpattern.image);
            I2 = pattern.image;
            C = imfuse(I1, I2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(C);
            hold on
            scatter(centers(1,:), centers(2,:),50, 'b.');
            scatter(centers(1,:)+qs(1,:), centers(2,:)+qs(2,:),50, 'r.');
            q = quiver(centers(1,:), centers(2,:), qs(1,:), qs(2,:), 'LineWidth', 1);
            q.AutoScale = 'off';
            hold off
        end
    end
end