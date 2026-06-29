function obj = betaClassic(obj, testPat, refPat, options)
%     if ~isempty(options.MakeGIF)
%         testPat.show_pattern;
%         xlabel(num2str(1))
%         frame = getframe(gcf);
%         im{1} = frame;
%     end
    if options.plotCC
        figure;
        obj.cctiles = tiledlayout(options.plotCC, 3,'TileSpacing','tight');
        rois = 1:obj.roi.numRois/options.plotCC:obj.roi.numRois;
    end
    rs = zeros(3, obj.roi.numRois);
    qs = zeros(3, obj.roi.numRois);
    testPC_p = obj.ft.pc2phosphorFrame(testPat.patternCenter);
    refPC_p = obj.ft.pc2phosphorFrame(refPat.patternCenter);
    rRef_p = obj.roiCenters_p - refPC_p;
    for i = 1:obj.roi.numRois
        ri_p = rRef_p(:,i);
        shiftRoi_p = obj.get_shift(testPC_p, refPC_p, ri_p, options.FGuess);
        shiftRoi = (obj.ft.phosphorFrame2imageVec(shiftRoi_p));
        [refRoi, rrange, crange] = obj.roi.get_squareRoi(...
            refPat, i ...
            );
        [testRoi, srrange, scrange] = ...
            obj.roi.get_squareRoiShifted( ...
                testPat, i, shiftRoi ...
        );
        ccMap = get_crossCorrelation(refRoi, testRoi, obj.roi, ...
            'Gradient',testPat.SimData.Gradient);
        q = obj.get_subpixShiftImage(real(ccMap));
        if options.plotCC && any(rois == i)
            obj.plotCC(real(ccMap), testRoi, refRoi, q)
        end
        q = q + [mean(scrange) - mean(crange);...
                 mean(srrange) - mean(rrange)];
        q_p = obj.ft.imageVec2phosphorFrame(q);
        qstar_p = obj.reverse_PC_offset(testPC_p, refPC_p, q_p, ri_p);
        qs(:,i) = qstar_p/norm(ri_p);
        rs(:,i) = ri_p/norm(ri_p);
    end
    obj.q = qs;
%     obj.useq = obj.getGoodIndicies(qs);
%     obj.calc_betaDetector(refPat, rs(:,obj.useq), qs(:,obj.useq));
%     obj.useq = logical(ones([1,size(qs,2)]));
    if any(mad(qs,1,2) > obj.MADThreshold)
        [~, outlierInd] = rmoutliers(qs');
        obj.useq = ~outlierInd';
    else
        obj.useq = true(1, size(qs,2)); %logical(ones([1,size(qs,2)]));
    end
    obj.calc_betaDetector(refPat, rs, qs);
    obj.g = obj.getg(refPat.g);
%     obj.calc_betaDetector(refPat, rs, qs);
%     if ~isempty(options.MakeGIF)
%         plot(obj, refPat)
%         xlabel(num2str(2))
%         frame = getframe(gcf);
%         im{2} = frame;
%         obj.makeGIF(im, options.MakeGIF)
%         close gcf
%     end
end