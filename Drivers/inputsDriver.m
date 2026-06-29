clc; clear; close all
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Define inputs
clc
% Filenames for inputs,
% If this cell array is not given a window will pop-up
% filenames = {'scanfile', 'C:\Users\wggilli\Documents\Al Void Nucleation\Reanalyzed Map Data 18.ctf'; ...
%              'imagefile', 'C:\Users\wggilli\Documents\Al Void Nucleation\ProcessedImages\0_0.tiff'};


% filenames = {'scanfile', '\\snl\Collaborative\hrebsdshare\Older Projects\DummySets\SmallSi\SmallSi Specimen 1 Site 1 Map Data 4.ctf'; ...
%              'imagefile', '\\snl\Collaborative\hrebsdshare\Older Projects\DummySets\SmallSi\SmallSi Specimen 1 Site 1 Map Data 4 Images\SmallSi Specimen 1 Site 1 Map Data 4_01.tiff'};

% 1. Setting analysis options
% 2. Create the inputs
%   - grain clean up/reconstruction
%   - read ebsd
%   - consolidate information
% 3. phase replacement (if needed)
% 4. Pattern center calibrations (Al only)
% 5. Real reference HREBSD
% 6. Applying hyrbid method (simulated HREBSD on reference patterns)
% 7. Post processing (stress/strain other reference frames)


% filenames = {'scanfile', '\\snl\Collaborative\hrebsdshare\2023\Noell\Al2219HREBSD_Reanalyzed.ctf';...
%             'imagefile', '\\snl\Collaborative\hrebsdshare\2023\Noell\Al2219HREBSD Specimen 1 Site 1 Map Data 2 Images\Al2219HREBSD Specimen 1 Site 1 Map Data 2_000001.tiff'};

% filenames = {'scanfile', "C:\Users\wggilli\Documents\ResolutionTests\Data\Deformed Specimen 1 Site 1 def_1_1024res.ctf";...
%             'imagefile', "C:\Users\wggilli\Documents\ResolutionTests\Data\Deformed Specimen 1 Site 1 def_1_1024res\ProcessedImages\0_0.tiff"};

filenames = {'scanfile', "D:\2024\Clones\G3\Step2\Step2 Specimen 1 Site 1 Map Data 1.ctf";...
            'imagefile', "D:\2024\Clones\G3\Step2\data\b39a3164-90fc-47a7-9133-93920a99ff31.ebsp"};
analysisOptions = {"numCores", 3, 'iterationLimit',7,...
    'numRois',48, 'simulationType', 'Dynamic',...
    'assumption', 'free-surface'};


% filenames = {'scanfile', 'D:\Dec_29_2022\far_screw.ang';...
%             'imagefile', 'D:\Dec_29_2022\far_screw.up2'};

% Analysis options 
% analysisOptions = {"numCores", 4, 'iterationLimit',20,...
%     'numRois',48, 'simulationType', 'Dynamic'};

% Grain options
grainOptions = {"misorientation", 5};

% Set up data object
input = mtexHREBSD('fileOptions', filenames,...
                    'analysisOptions', analysisOptions,...
                    'grainOptions', grainOptions);
%%
% input = mtexHREBSD('analysisOptions', analysisOptions,...
%                     'grainOptions', grainOptions);

% This line replaces the Westgreen version of FeCu2Al7
input.replace_phase('FeCu2Al7 (Westgreen)', input.ebsd.CSList{5});

% ni_cs = crystalSymmetry('cubic',[3.523,3.523,3.523],'mineral','Nickel');
% for phase = input.ebsd.CSList
%     if class(phase{1}) == "crystalSymmetry"
%         input.replace_phase(phase{1}.mineral, ni_cs);
%     end
% end


% clear filenames analysisOptions grainOptions
%% Pattern center calibration
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

% pcs = zeros(length(indiciesPCCal), 3);
% for i = 1:length(indiciesPCCal)
%     testPat = input.get_pattern(indiciesPCCal(i));
%     pcs(i,:) = testPat.patternCenter;
% end
% figure
% scatter3(pcs(:,1), pcs(:,3), pcs(:,2))
% xlabel('x*')
% ylabel('z*')
% zlabel('y*')
%% Troubleshooting
clc 
C = parallel.pool.Constant(input);
parfor i = 1
    disp(class(C.Value.patterns))
end
%% Perform PCCal
% clc
% PCCal = mtexHREBSD_patternCenterCalibration(input,'indicies',indiciesPCCal);
% PCCal = PCCal.get_patternCenterCal(input);
% PCCal.plot(PCCal.get_patternCenterOffset);
% clc
input.patternCenterOffset = PCCal.get_patternCenterOffset;
% save('C:\Users\wggilli\Documents\Al Void Nucleation\PCCalInput.mat', "input")
% save('C:\Users\wggilli\Documents\Al Void Nucleation\PCCal.mat', "PCCal")

%%
load('C:\Users\wggilli\Documents\Al Void Nucleation\PCCal.mat')

%%
clc
figure
subset = get_subset(9,5,6,6,input.ebsd);
plot(input.ebsd(subset))
%%
clc
% output = hrebsdMainSimulated(input, 'Subset', subset);
% save('C:\Users\wggilli\Documents\Al Void Nucleation\PCCalKinematicSimulatedOutput.mat', "output")

input.analysis.simulationType = 'Dynamic';
input.analysis.Gradient = 0;
% output = hrebsdMainClassic(input, 'Subset', subset);
output = hrebsdMainSimulated(input);
% save('C:\Users\wggilli\Documents\Al Void Nucleation\PCCalOutputUpdated.mat', "output")
%%
load('C:\Users\wggilli\Documents\Al Void Nucleation\PCCalDynamicSimulatedOutput.mat', "output")
%%
updatedOutput = mtexHREBSD_hybrid(input, output, 1, 10);
%%
post = mtexHREBSD_postProcessing(input, output)
%%
post.multiplotStrainBeta(post,...
    'climsStrain', [1,1]*-0.0016 + [-1, 1]*0.0107,...
    'climsBeta', [-1,1]*5e-3 ...
    );
%%
post.get_multiplot('strain','refFrame', 'sample', 'clims', [-1,1]*5E-3);
%%
ebsdSubset = input.ebsd(subset);
subset2 = ebsdSubset('Aluminium').id;
strain = post.strain.phosphor(subset2);
% mean(strain)
mode(strain.M(1,3,:))
%%
clc
testingGUI(input)
%% Run HREBSD (Real reference patterns)
clc
% output = hrebsdMainClassic(input);

%% Run HREBSD (Simulated reference patterns)
clc
subset = get_subset(5,9,6,6,input.ebsd);
plot(input.ebsd(subset))
%%
% output = hrebsdMainSimulated(input, 'Subset', subset);
output2 = hrebsdMainClassic(input, 'Subset', subset);
save('C:\Users\wggilli\Documents\Al Void Nucleation\PCCalClassicOutput.mat', "output2")
%%
post = mtexHREBSD_postProcessing(input, output);
postSubset = post.do_subset(subset);
%%
idsCuAl2 = postSubset.hrebsd.ebsd.id(postSubset.hrebsd.ebsd.grainId == 111);
idsAl = postSubset.hrebsd.ebsd.id(postSubset.hrebsd.ebsd.grainId ~= 111);
postSubset.multiplotStrainBeta(postSubset, ...
    'clims', [-1,1]*5e-3,...
    'refFrame', 'sample', ...
    'noMaskIds', idsCuAl2)
postSubset.multiplotStrainBeta(postSubset, ...
    'clims', [-1,1]*5e-2,...
    'refFrame', 'sample', ...
    'noMaskIds', idsAl)
%%
postSubset.get_multiplot('beta', 'doScaleBar', 0, 'refFrame', 'sample', 'clims',[-1,1]*3E-3, 'FigSize',[15,15])
%% F-Delta
output = hrebsdMainFDelta(input);
%%
ind = 38381;
grainId = input.ebsd.grainId(ind);
refId = input.refIds(grainId);
testPat = input.get_pattern(ind);
refPat = input.get_pattern(refId);
simPat = testPat.get_simulatedPattern(input);
% refPat.show_pattern
testPat.show_pattern
simPat.show_pattern
%%
close all; clc
ind = 22225; % 32367; %38381;
testPat = input2.get_pattern(ind);
testPat.show_pattern
test = FDeltaCalcF(testPat, input, 'TroubleshootingPlots', 1, 'ConvergencePlots', 1);
test2 = dynamicSimulatedCalcF(testPat, input, 'TroubleshootingPlots', 1, 'ConvergencePlots', 1);

%%
test = get_subset(5,9,6,5,input.ebsd);
plot(input.ebsd(test))

function subsetIds = get_subset(xmin, ymin, width, height, ebsd)
    x = ebsd.prop.x;
    y = ebsd.prop.y;
    xLogic = x >= xmin & x <= xmin + width;
    yLogic = y >= ymin & y <= ymin + height;
    subsetLogic = find(xLogic.*yLogic);
    subsetIds = ebsd.id(subsetLogic);
end