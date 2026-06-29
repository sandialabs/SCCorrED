tensorPath = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\data\tensor';
EMdataPath = '\Users\wggilli\Documents\EMsoft-5.0.20230206.-Win64\data';
EMsoftPath = '\Users\wggilli\Documents\EMsoft-5.0.20230206.-Win64\bin\';
materialPath = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\data\Materials';

path = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\HelperFunctions';
fname = fullfile(path, 'mtexHREBSD_paths');
save(fname, 'tensorPath', 'EMdataPath', 'EMsoftPath', 'materialPath');

%%
configFile = fullfile( ...
    getenv("SYSTEMDRIVE"), ...
    "Users", ...
    getenv('username'), ...
    ".config\mtexHREBSD\mtexHREBSDConfig.json" ...
    );
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
val = jsondecode(str);
