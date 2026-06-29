function [F, beta] = get_grainF(mtexHREBSD)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
numPoints = length(mtexHREBSD.ebsd);
beta = zeros(3,3,numPoints);
F = zeros(3,3,numPoints);
viableRefs = viableRefIds(mtexHREBSD,20);

end