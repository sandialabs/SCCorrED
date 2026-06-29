function ccMap = get_crossCorrelation(ref_roi, test_roi, roiSettings, varargin)
%get_crossCorrelation Summary of this function goes here
%   Detailed explanation goes here
    args = get_args(varargin{:});
%     args
    if args.Gradient
        ref_roi_filtered = roiSettings.windowfunc.*ref_roi;
        test_roi_filtered = roiSettings.windowfunc.*test_roi;
        [ref_Gx, ref_Gy] = imgradientxy(ref_roi_filtered);
        [test_Gx, test_Gy] = imgradientxy(test_roi_filtered);
        ref_G = complex(ref_Gx, ref_Gy);
        test_G = complex(test_Gx, test_Gy);
        ref_fft = fftn(ref_G);
        test_fft = fftn(test_G);
    else
        ref_fft = fftn(roiSettings.windowfunc.*ref_roi);
        test_fft = fftn(roiSettings.windowfunc.*test_roi);
    end
    ccMap = fftshift(ifftn(roiSettings.custfilt.*ref_fft.*conj(roiSettings.custfilt.*test_fft)));
%     ccMap = fftshift(ifftn(roiSettings.custfilt.*test_fft.*conj(roiSettings.custfilt.*ref_fft)));
    if args.Plot
        if args.Gradient
            plot_CC_gradient(ref_G, test_G, ccMap)
        else
            plot_CC(ref_roi, test_roi, ccMap)
        end
    end
end


function p = create_inputParser
    p = inputParser;
    addParameter(p, 'Gradient', 0)
    addParameter(p, 'Plot', 0)
end


function args = get_args(varargin)
    p = create_inputParser;
    parse(p, varargin{:});
    args = p.Results;
end


function plot_CC(ref_roi, test_roi, CC)
    figure
    tiledlayout(1,3, "TileSpacing","compact", "Padding","compact")
    ax = nexttile;
    imagesc(ref_roi); axis image; axis off; colormap(ax, "gray");
    ax = nexttile;
    imagesc(test_roi); axis image; axis off; colormap(ax, "gray");
    ax = nexttile;
    imagesc(CC); axis image; axis off; colormap(ax, "jet");
end


function plot_CC_gradient(refG, testG, CC)
    figure
    tiledlayout(2,3, "TileSpacing","compact", "Padding","compact")
    ax = nexttile;
    imagesc(real(refG)); axis image; axis off; colormap(ax, "gray");
    ax = nexttile;
    imagesc(real(testG)); axis image; axis off; colormap(ax, "gray");
    ax = nexttile([2,1]);
    imagesc(real(CC)); axis image; axis off; colormap(ax, "jet");

    ax = nexttile;
    imagesc(imag(refG)); axis image; axis off; colormap(ax, "gray");
    ax = nexttile;
    imagesc(imag(testG)); axis image; axis off; colormap(ax, "gray");
%     ax = nexttile;
%     imagesc(imag(CC)); axis image; axis off; colormap(ax, "jet");

%     figure
%     imagesc(real(CC)); axis image; axis off; colormap jet 
end