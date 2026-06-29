function ccMap = get_xcorr2(ref_roi, test_roi, roiSettings, options)
%get_crossCorrelation Summary of this function goes here
%   Detailed explanation goes here
%     args = get_args(varargin{:});
%     args
    ref_filt = ...
        ifftn(roiSettings.custfilt.*fftn(roiSettings.windowfunc.*ref_roi));
    test_filt = ...
        ifftn(roiSettings.custfilt.*fftn(roiSettings.windowfunc.*test_roi));
    if options.NormXCorr
        ccMap = abs(normxcorr2(ref_filt, test_filt));
    else
        ccMap = xcorr2(ref_filt, test_filt);
    end
end