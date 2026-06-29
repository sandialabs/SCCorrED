function ids = singlePixelGrains(ebsd, grains)
    ids = zeros(size(grains.id));
    for i = 1:size(grains.id)
        if sum(ebsd.grainId == i) <= 1
            ids(i) = 1;
        end
    end
end