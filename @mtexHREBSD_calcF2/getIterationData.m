function getIterationData(obj, pat, iter, options)
    if options.SaveIterations
        obj.ImageIterations(:,:,iter) = uint8(rescale(pat.image,0,255));
        obj.qIteration(:,:,iter) = obj.q;
    end
end