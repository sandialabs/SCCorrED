%=========================================================================%
% This script runs mtexHREBSD using real reference patterns and saves the 
% results as .mat files for separate processing.
%
% By Thomas Bennett, 2025-04-22
%
% This script uses MTEX and mtexHREBSD.
%
% Updates
%  2025-07-23: - refined comments and formatting
%              - added example post-processing and plotting
%
%=========================================================================%

%% Setup
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% Get date and time
dt = datetime("now");
dtStr = string(datetime(dt,"Format","yyyyMMdd_HHmm"));

%% Define inputs

suffix = "real_ref"; % suffix for output file names

% specify files to import
filenames = {'scanfile', "<PATH TO ORIENTATION DATA (e.g., .ang, .ctf, .h5oina)>";...
             'imagefile', "<PATH TO PATTERN FILE (.up1/2, .h5oina) OR FIRST .tif IMAGE>"};

% specify HREBSD analysis options
analysisOptions = {"numCores", 4, 'iterationLimit',15,...
    'numRois',48, 'simulationType', 'Dynamic',...
    'assumptions', "trace=0"};
% note: for real reference patterns, the options 'iterationLimit' and 
%       'simulationType' may be set arbitrarily

% specify grain segmentation options
grainOptions = {"misorientation", 5, "minSize", 5};

% mtexHREBSD input object
input = mtexHREBSD('fileOptions', filenames,...
                    'analysisOptions', analysisOptions,...
                    'grainOptions', grainOptions);

% adjust reference points manually if necessary
% input.refIds = ...;

%% Export reference patterns

% make a folder to save reference patterns
pat_dir = ".\reference_patterns";
if ~isfolder(pat_dir);  mkdir(pat_dir);  end

% do not save reference patterns for grains smaller than this minimum
minGrainSize = 5; % pixels

% loop through reference points and save patterns
for i=1:length(input.refIds)
    grainSize = input.grains(i).grainSize;
    % only save references for sufficiently large grains
    if grainSize >= minGrainSize
        % import the pattern
        pat = input.get_pattern(input.refIds(i));
        % rescale and save the pattern
        imwrite(rescale(pat.image),fullfile(pat_dir, ...
                sprintf("ref_pat_%d.png",input.refIds(i))));
    end
end

%% Examine reference points

fig = figure;
% plot band contrast
plot(input.ebsd,input.ebsd.prop.bc);
mtexColorMap black2white;
hold on;
% plot grain boundaries
plot(input.grains.boundary);
% plot reference points
plot(input.ebsd.prop.x(input.refIds),input.ebsd.prop.y(input.refIds), ...
    'ro','MarkerFaceColor','r','MarkerSize',7);
exportgraphics(fig,strcat("ref_points_map_",suffix,".png"), ...
    "Resolution", 200);

%% Check calculations on a single reference pattern

index = 1; % user specified index
refId = input.refIds(input.ebsd(index).grainId);

% generate a figure to show the test and reference points
figure;
plot(input.ebsd);
hold on;
plot(input.ebsd(index).prop.x, input.ebsd(index).prop.y, 'r.'); % test
plot(input.ebsd(refId).prop.x, input.ebsd(refId).prop.y, 'b.'); % ref

testPat = input.get_pattern(index);
refPat = input.get_pattern(refId);
figure; testPat.show_pattern;
figure; refPat.show_pattern;

% calculate lattice distortion
calcfstrain = input.initializeCalcF('AnalysisType', 0);
solution = calcfstrain.getBeta(testPat,refPat);
% plot the reference pattern with measured shifts
figure;
plot(calcfstrain,refPat,'multiplier',5,'showFit',1);
beta = solution.beta % show the elements of beta in the phosphor frame

%% HREBSD

% initialize a calcF object
calcfstrain = input.initializeCalcF('AnalysisType', 0);

% log HREBSD settings
input.iniWrite("File",strcat(".\HREBSD_settings_",dtStr,".txt"));

% run HREBSD processing
tic
output = mtexHREBSDMain(input,...
    'calcFStrain', calcfstrain, ...
    'ParallelScheme', 2, ...
    'numCores', 4);
executionTime = toc;

% save data to .mat files
save(strcat("results_",suffix),"output","input");
save(strcat("workspace_",suffix));

%% Post-processing

% load in processed HREBSD data (optional)
% load("<PATH TO 'results_.mat' FILE EXPORTED ABOVE>","input","output");

% run post processing
post = mtexHREBSD_postProcessing(input, output);

% set colormap limits
rangeRot = 1e-1;
rangeStrain = 2e-2;
rangeStress = 500e6; % Pascals

% plot lattice distortion (beta) in the phosphor reference frame
post.get_multiplot("beta", "refFrame", "phosphor", "doSmooth", 0, ...
                   "clims", [-1,1]*rangeRot);

% plot strain in the sample frame
post.multiplotStrainBeta(post,'climsStrain', [-1,1]*rangeStrain, ...
                         'climsBeta', [-1,1]*rangeRot,'refFrame', ...
                         'sample','doSmooth', 0);

% plot stress in the sample frame
post.get_multiplot("stress", "refFrame", "sample", "doSmooth", 0, ...
                   "clims", [-1,1]*rangeStress);

% plot sum of squares error (SSE)
fig = figure;
plot(post.hrebsd.ebsd,post.data.SSE);
mtexColorbar
mtexColorMap WhiteJet01
clim([0,0.05]);