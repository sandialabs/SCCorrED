function refId = grainRefId(ebsdGrain, grain)
    dists = distFromGB(ebsdGrain, grain);
    distsNormalized = normalize0to1(dists);
    GROD = angle(ebsdGrain.orientations, grain.meanOrientation);
    GRODNormalized = normalize0to1(1./GROD);
    if isfield(ebsdGrain.prop, 'bc')
        qualNormalized = normalize0to1(ebsdGrain.bc);
    else
        qualNormalized = normalize0to1(ebsdGrain.iq);
    end
    metric = 0.3333.*distsNormalized +...
        0.3333.*GRODNormalized + 0.3333*qualNormalized;
    [~,ind] = max(metric);
    refId = ebsdGrain.id(ind);
end


function xnormalized = normalize0to1(x)
    xnormalized = (x - min(x))./(max(x) - min(x));
end