clc; clear; close all;
savedir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting';
load(fullfile(savedir,"test_angStep3_2.mat"));

%%
plot(test_angStep3_2.grains)
roi = drawrectangle
xRange = [roi.Position(1), roi.Position(1)+roi.Position(3)];
yRange = [roi.Position(2), roi.Position(1)+roi.Position(4)];

%%
ebsd = test_angStep3_2.ebsd;
subScanInds = zeros(length(ebsd),1);
xLogic = xRange(1) <= ebsd.prop.x & ebsd.prop.x <= xRange(2);
yLogic = yRange(1) <= ebsd.prop.y & ebsd.prop.y <= yRange(2);
logic = xLogic & yLogic;
subScan = ebsd(logic);
plot(test_angStep3_2.grains)
hold on
plot(subScan)
%%
clc
subScanIds = ebsd.id(logic);
data = cell(length(subScanIds),1);
parfor i = 1:length(subScanIds)
    temp = test_angStep3_2;
    id = subScanIds(i);
    testPat = temp.get_pattern(id);
    try
        data{i} = dynamicSimulatedCalcF(testPat, temp);
    catch ME
        disp(ME)
    end
    disp(i)
end
%%
test_angStep3_2 = test_angStep3PCCal;
subsetAnalysis = mtexHREBSD_main(test_angStep3_2, 1);
subscanIds = zeros(length(data),1);
for i = 1:length(data)
    try
        scanIndex = data{i}.scanIndex;
        subscanIds(i) = scanIndex;
        subsetAnalysis.g(:,:,scanIndex) = data{i}.g;
        subsetAnalysis.beta(:,:,scanIndex) = data{i}.beta;
        subsetAnalysis.F(:,:,scanIndex) = data{i}.F;
    catch me
        disp(i)
        subscanIds(i) = scanIndex + 1;
    end
end
%%
postDyn = mtexHREBSD_postProcessing(test_angStep3_2, subsetAnalysis);
%%
postDyn = postDyn.do_subset(subscanIds);
%%
postDyn.get_multiplot('beta', 'clims', [-1,1]*2.5E-3, 'doSmooth', 0, ...
    'refFrame','phosphor')
postDyn.get_multiplot('strain', 'clims', [-1,1]*2.5E-3, 'doSmooth', 0, ...
    'refFrame','phosphor')
postDyn.get_multiplot('stress', 'clims', [-1,1]*5E8, 'doSmooth', 0, ...
    'refFrame','sample')
%%
grainId = 114;
grainStep3 = test_angStep3_2.get_grain(grainId);
grainTest = mtexHREBSD_main(grainStep3);