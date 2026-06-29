function distFromCentroid = distFromBoundary(obj, grainId)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     for i = 1:length(a.grains)
       grain_i = obj.grains(grainId);
       ebsd_i = obj.ebsd(obj.ebsd.grainId == grain_i.id);
       centroid_i = grain_i.centroid;
       distFromCentroid = zeros(length(ebsd_i), 1);
       for j = 1:length(ebsd_i)
           deltaX_j = ebsd_i.prop.x(j) - centroid_i(1);
           deltaY_j = ebsd_i.prop.y(j) - centroid_i(2);
           distFromCentroid(j) = mean((deltaX_j.^2 - deltaY_j.^2).^0.5);
       end
%     end
%     figure
%     plot(ebsd_i);
%     hold on 
%     scatter(obj.ebsd.prop.x(obj.refIds(grain_i.id)), obj.ebsd.prop.y(obj.refIds(grain_i.id)))
%     scatter(obj.ebsd.prop.x(minId), obj.ebsd.prop.y(minId), 'kx')
%     scatter(grain_i.boundary.x, grain_i.boundary.y)
%     out = 1;
end