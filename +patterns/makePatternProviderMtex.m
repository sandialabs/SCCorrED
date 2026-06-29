function patternProvider = makePatternProviderMtex(ebsd, first_image)

[~, ~, ext] = fileparts(first_image);

    switch ext
        case '.h5'
            patternProvider =...
                patterns.H5PatternProvider(first_image);
        case '.h5oina'
            patternProvider =...
                patterns.H5oinaPatternProvider(first_image);
        case {'.up1', '.up2'}
            patternProvider =...
                patterns.UPPatternProvider(first_image);
        case { '.jpg', '.jpeg', '.tif', '.tiff', '.bmp', '.png'}
            x = single(ebsd.prop.x);
            y = single(ebsd.prop.y);
            x(abs(x) < 1E-5) = 0;
            y(abs(y) < 1E-5) = 0;
            X = unique(x);
            Y = unique(y);
            xStep = X(2) - X(1);
            yStep = Y(2) - Y(1);
            Nx = length(X);
            Ny = length(Y);
            patternProvider =...
                patterns.ImagepatternProvider(...
                first_image,...
                'Square',...
                length(ebsd),...
                [Nx, Ny],...
                [ebsd.prop.x(1), ebsd.prop.x(2)],...
                [xStep, yStep]);
            
        otherwise
            error('OpenXY:PatternProvider:UnrecognizedExtention', ...
                'Unrecognized file extention %s', ext)
    end
end

