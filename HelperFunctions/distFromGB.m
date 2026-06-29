function dists = distFromGB(ebsd, grain)
    arguments
        ebsd EBSD
        grain grain2d
    end
    xebsd = repmat(ebsd.prop.x, 1, size(grain.x,1));
    yebsd = repmat(ebsd.prop.y, 1, size(grain.y,1));
    xgrain = repmat(grain.x', size(ebsd.prop.x,1), 1);
    ygrain = repmat(grain.y', size(ebsd.prop.y,1), 1);
    dx = xebsd - xgrain;
    dy = yebsd - ygrain;
    alldists = (dx.^2 + dy.^2).^0.5;
    dists = min(alldists,[],2);
end
