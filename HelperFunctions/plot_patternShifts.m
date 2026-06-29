function plot_patternShifts(pattern, shifts, roi)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    image = uint8(rescale(pattern.image, 0, 255));
    centers = roi.centers;
    f = figure;
    imshow(image);
    hold on
    quiver(centers(1,:), centers(2,:), shifts(1,:), shifts(2,:), 'rx');
end