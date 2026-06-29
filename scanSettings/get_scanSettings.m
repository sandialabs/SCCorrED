function scanSettings = get_scanSettings(varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    if ~isempty(varargin{:})
        options = varargin{:};
    else 
        options = {};
    end
    check_for_scan = check_scanfile(options);
    current_row = size(options, 1);
    if ~check_for_scan
        current_row = current_row + 1;
        options{current_row, 1} = 'scanfile';
        scanfile = get_scanfilename;
        options{current_row, 2} = scanfile;
    end
    check_for_image = check_imagefile(options);
    if ~check_for_image
        current_row = current_row + 1;
        options{current_row, 1} = 'imagefile';
        options{current_row, 2} = get_imagefilename;
    end
    options_struct = convert_to_struct(options);
    [~,~,ext] = fileparts(options_struct.scanfile);
    switch ext
        case '.ang'
            scanSettings = mtexHREBSD_scanSettings_ang(options_struct);
        case '.ctf'
            scanSettings = mtexHREBSD_scanSettings_ctf(options_struct);
        case '.h5oina'
            scanSettings = mtexHREBSD_scanSettings_h5oina(options_struct);
    end
end


function options_struct = convert_to_struct(options)
    options_struct = struct();
    for i = 1:size(options, 1)
        options_struct.(options{i, 1}) = options{i, 2};
    end
end

function check = check_imagefile(options)
    ind = find(strcmp(options, 'imagefile'));
    if ~isempty(ind)
        check = ind;
    else 
        check = 0;
    end
end


function imagefile = get_imagefilename
    [name, path, ~] = uigetfile({
        '*.jpg;*.jpeg;*.tif;*.tiff;*.bmp;*.png','Image Files'
        '*.up1;*.up2', 'OIM Uncompressed Pattern Format'
        '*.ebsp','EBSP Format'
        },...
        'Select the First Image of the Scan or pattern archive');
    imagefile = fullfile(path, name);
end


function check = check_scanfile(options)
    ind = find(strcmp(options, 'scanfile'));
    if ~isempty(ind)
        check = ind;
    else 
        check = 0;
    end
end


function scanfile = get_scanfilename
    [name, path, ~] = uigetfile({
        '*.ang;*.ctf', 'Scan Files (*.ang,*.ctf)'
        '*.h5', 'OIM HDF5 Files (*.h5)'
        },'Select a Scan File');
    scanfile = fullfile(path, name);
end

