function subsetIds = get_subset(xmin, ymin, width, height, ebsd)
    x = ebsd.prop.x;
    y = ebsd.prop.y;
    xLogic = x >= xmin & x <= xmin + width;
    yLogic = y >= ymin & y <= ymin + height;
    subsetLogic = find(xLogic.*yLogic);
    subsetIds = ebsd.id(subsetLogic);
end