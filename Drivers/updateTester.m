clc; clear; close all;
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
%%
% clear input; clc
ANALYSISOPTS = {"numCores", 1, 'iterationLimit',8,...
    'numRois',48, 'simulationType', 'Dynamic',...
    'assumption', 'free-surface'};
GRAINOPTS = {"misorientation", 5};
% 
DATADIR = "C:\Users\wggilli\Documents\ResolutionTests\Data";
FILENAME = "Deformed Specimen 1 Site 1 def_1_1024res.h5oina";

% DATADIR = "\\snl\Collaborative\hrebsdshare\2023\DeitzHiperco";
% FILENAME = "APDCheck Specimen 1 Site 6 Map Data 4.ctf";

input2 = createInput(ANALYSISOPTS, GRAINOPTS, fullfile(DATADIR, FILENAME));
% input.patternCenterOffset = os;
% 
% DATADIR = "D:\Dec_29_2022";
% FILENAME = "far_screw.ang";
% input = createInput(ANALYSISOPTS, GRAINOPTS, fullfile(DATADIR, FILENAME));
%%
clc
% Custom PCCal points using only one phase using get_grid function
spacing = 2;
% phase = 'Aluminium';
phase = 'Ni';
indiciesPCCal = get_grid(input.ebsd, spacing, 'phase', phase);
figure
plot(input.ebsd);
hold on

scatter(input.ebsd.prop.x(indiciesPCCal),input.ebsd.prop.y(indiciesPCCal),'rx');
legend('location','bestoutside')
%%
PCCal = mtexHREBSD_patternCenterCalibration(input,'indicies',indiciesPCCal);
PCCal = PCCal.get_patternCenterCal(input);
%%
path = "C:\Users\wggilli\Documents\hrebsdTestData";
fname = fullfile(path, "inPCCal.mat");
if isfile(fname)
    load(fname)
else
    input.patternCenterOffset = PCCal.get_patternCenterOffset; 
    save("C:\Users\wggilli\Documents\hrebsdTestData\inPCCalh5.mat", "input");
end

%% Updated main testing
clc
fname = fullfile(path, "outRealh5.mat");
if isfile(fname)
    load(fname)
else
    clear calcfstrain calcfalpha test
    calcfstrain = input.initializeCalcF('AnalysisType',0);
    calcfalpha = input.initializeCalcF;
    outReal = mtexHREBSDMain(input, ...
        'calcFStrain', calcfstrain,...
        'calcFAlpha',calcfalpha,...
        'ParallelScheme',2,...
        'numCores', 5,...
        'ChunkSize',10, ...
        'Verbose',1);
    save( ...
        "C:\Users\wggilli\Documents\hrebsdTestData\outRealh5.mat", "outReal" ...
        );
end
%%
fname = fullfile(path, "outSimh5.mat");
if isfile(fname)
    load(fname)
else
    calcfstrain = input.initializeCalcF('AnalysisType',1);
    outSim = mtexHREBSDMain(input, ...
        'calcFStrain', calcfstrain,...
        'ParallelScheme',2,...
        'numCores', 5,...
        'ChunkSize',10, ...
        'Verbose', 1);
    save("C:\Users\wggilli\Documents\hrebsdTestData\outSimh5.mat", "outSim");
end
%%
clear calcfstrain calcfstrain2
calcfstrain = input.initializeCalcF('AnalysisType',2);
calcfstrain2 = input2.initializeCalcF('AnalysisType',2);
p = input.get_pattern(1);
p2 = input2.get_pattern(1);
p2.patternCenter = p.patternCenter;
% r = input.get_pattern(subset);
calcfstrain.getBeta(p, 'MakeGIF', "C:\Users\wggilli\Documents\hrebsdTestData\testctf.gif")
calcfstrain2.getBeta(p2, 'MakeGIF', "C:\Users\wggilli\Documents\hrebsdTestData\testh5.gif")
%%
% clc
% subset = input.refIds;
calcfstrain = input.initializeCalcF('AnalysisType',2)
calcfstrain2 = input2.initializeCalcF('AnalysisType',2)
p = input.get_pattern(1);
p2 = input2.get_pattern(1);
p2.patternCenter = p.patternCenter;
% r = input.get_pattern(subset);
calcfstrain.getBeta(p2, 'ConvergencePlots', 1, 'MakeGIF', "C:\Users\wggilli\Documents\hrebsdTestData\test.gif")
% calcfstrain.AnalysisType = 0;
% calcfstrain.getBeta(p,r,'MakeGIF', "C:\Users\wggilli\Documents\hrebsdTestData\classic.gif")
%%
tester = input.initializeCalcF('AnalysisType', 1);
p = input.get_pattern(randi(length(input.ebsd), 1))
tester.getBeta(p, ...
    'ConvergencePlots', 1, ...
    'MakeGIF', "C:\Users\wggilli\Documents\hrebsdTestData\test.gif" ...
    );
%%
calcfstrain = input.initializeCalcF('AnalysisType',1);
outRefs = mtexHREBSDMain(input, ...
    'calcFStrain', calcfstrain,...
    'ChunkSize',1, ...
    'Subset', subset,...
    'Verbose', 1, ...
    'Subset', subset);
%%
calcfstrain = input.initializeCalcF('AnalysisType',0);
test2 = mtexHREBSDMain(input, ...
    'calcFStrain', calcfstrain,...
    'ParallelScheme',2,...
    'numCores', 4,...
    'ChunkSize',20, ...
    'Verbose',1);
%%
clc
F(1) = parfeval(@bingus,1);
F(2) = parfeval(@bongus,1);
try
    [i,o,k,r] = readFutures(F(1), 1, 1);
catch me
    if strcmp(me.identifier, ...
            'MATLAB:parallel:future:FetchNextFutureErrored')
        disp('Do the thing')
    end
end

function [i, out, k, range] = readFutures(F, k, chunk)
    [i, out] = fetchNext(F);
    k = k + 1;
    lo = (k-1)*chunk + 1;
    hi = (k-1)*chunk + size(out,1);
    range = lo:hi;
    nUpdateWaitbar;
end


function out = bingus
    out = 0;
    error("This is the error")
end
function out = bongus
    out = 1;
end

function hrebsdinput = createInput(analysisopts, grainopts, file)
    [path,fname,ext] = fileparts(file);
    if strcmp(ext, '.ctf')
        scanpath = fullfile(path, fname+ext);
        imgpath = fullfile(path, fname,"ProcessedImages","0_0.tiff");
%         imgpath = fullfile(path, fname+" Images",fname+"_00001.tiff");
    elseif strcmp(ext, '.h5oina')
        scanpath = file;
        imgpath = file;
    else
        scanpath = fullfile(path, fname+ext);
        imgpath = fullfile(path, fname+".up2");
    end
    filenames = {'scanfile', scanpath; 'imagefile', imgpath};
    hrebsdinput = mtexHREBSD('fileOptions', filenames,...
                    'analysisOptions', analysisopts,...
                    'grainOptions', grainopts);
end
