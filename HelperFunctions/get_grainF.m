function [F, beta] = get_grainF(mtexHREBSD)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
numPoints = length(mtexHREBSD.ebsd);
beta = zeros(3,3,numPoints);
F = zeros(3,3,numPoints);
viableRefs = viableRefIds(mtexHREBSD.ebsd, 10);
refPattern = mtexHREBSD.get_pattern(viableRefs.firstLocal);
for i = 1:numPoints
    testPattern = mtexHREBSD.get_pattern(i);
    calcF_i = classicCalcF(refPattern, testPattern, mtexHREBSD);
    beta(:,:,i) = calcF_i.beta;
    F(:,:,i) = calcF_i.F;
end