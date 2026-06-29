clc; clear; close all;

savedir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting';
dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Away_screw";
hrebsdAnalysisFilename = "HREBSD_gUpdate2";
minGrainSize = 3;
hrebsdSteps = [0,1,2,3];
dicSteps = [2,3,4,5];

grainsToUseFDelta = {[146,157,142,116,110,118,117,108,139,149,150,129,114,106,145,98,107,112,94], ...
               [146,161,143,117,109,121,120,116,104,144,153,130,115,107,112,108,95], ...
               [151, 164, 148, 122, 113, 124, 123,108,120,149,156,135,118,111,117,112,99], ...
               143};%[113,101,141,143,156,149,127,112,140,114,105,110,93,106,117,115]};

for i = 1:length(hrebsdSteps)
    registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg";
    hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal";
    dicFilename = "away_screw_00"+num2str(dicSteps(i));
    reg = load(fullfile(savedir, registrationFilename)).(registrationFilename);
    hrebsdInput = load(fullfile(savedir, hrebsdFilename)).(hrebsdFilename);
    dic = load(fullfile(dicdir, dicFilename));
%     subsetIndicies = findDataSubset(reg,dic,hrebsdInput);
    subsetIndicies = findDataSubsetFullGrains(reg,dic,hrebsdInput,minGrainSize);
%     subsetIndicies = [];
%     grains = grainsToUseFDelta{i};
%     for j = 1:length(grains)
%         subsetIndicies = [subsetIndicies; hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == grains(j))];
%     end
    figure
    plot(hrebsdInput.ebsd(subsetIndicies))
    hrebsdInput.analysis.iterationLimit = 10;
    hrebsdOutput = hrebsdMainClassic(hrebsdInput, 'Subset', subsetIndicies);
    saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_classic.mat");
    save(saveTo, "hrebsdOutput");
    grainsToUse = unique(hrebsdInput.ebsd.grainId(subsetIndicies));
    calcFRefs = cell(length(grainsToUse), 1);
    h = waitbar(0, "Starting reference patterns", 'Name', "Step "+num2str(hrebsdSteps(i)));
    nUpdateWaitbar(length(calcFRefs), h);
    C =  parallel.pool.Constant(hrebsdInput);
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    parfor j = 1:length(calcFRefs)
        try
            refPat = C.Value.get_pattern(C.Value.refIds(grainsToUse(j)));
            calcFRefs{j} = {grainsToUse(j),dynamicSimulatedCalcF(refPat, C.Value)};
        catch me
            disp(j)
            disp(me)
        end
        send(D, 1)
    end
    saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_calcFRefs.mat");
    save(saveTo, "calcFRefs");
    close(h)
end
% i = 4;
% subsetIndicies = [];
% grains = grainsToUseFDelta{i};
% for j = 1:length(grains)
%     subsetIndicies = [subsetIndicies; hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == grains(j))];
% end
% hrebsdOutputFDelta = hrebsdMainFDelta(hrebsdInput, 'Subset', subsetIndicies);
% saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_FDelta.mat");
% save(saveTo, "hrebsdOutputFDelta");
%% From save data
clc; clear; close all;

savedir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting';
dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Away_screw";
hrebsdAnalysisFilename = "HREBSD_gUpdate2";
minGrainSize = 50;
hrebsdSteps = [0,1,2,3];
dicSteps = [2,3,4,5];


for i = 1:length(hrebsdSteps)
    registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg";
    hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal";
    hrebsdOutputFilename = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_classic.mat");
    dicFilename = "away_screw_00"+num2str(dicSteps(i));
    reg = load(fullfile(savedir, registrationFilename)).(registrationFilename);
    hrebsdInput = load(fullfile(savedir, hrebsdFilename)).(hrebsdFilename);
    dic = load(fullfile(dicdir, dicFilename));
%     subsetIndicies = findDataSubset(reg,dic,hrebsdInput);
    subsetIndicies = findDataSubsetFullGrains(reg,dic,hrebsdInput,minGrainSize);
%     subsetIndicies = [];
%     grains = grainsToUseFDelta{i};
%     for j = 1:length(grains)
%         subsetIndicies = [subsetIndicies; hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == grains(j))];
%     end
    hrebsdOutput = load(hrebsdOutputFilename).("hrebsdOutput");
    figure
    plot(hrebsdInput.ebsd(subsetIndicies))
    drawnow

    grainsToUse = unique(hrebsdInput.ebsd.grainId(subsetIndicies));
    calcFRefs = cell(length(grainsToUse), 1);
    h = waitbar(0, "Starting reference patterns", 'Name', "Step "+num2str(hrebsdSteps(i)));
    nUpdateWaitbar(length(calcFRefs), h);
    C =  parallel.pool.Constant(hrebsdInput);
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    parfor j = 1:length(calcFRefs)
        try
            refPat = C.Value.get_pattern(C.Value.refIds(grainsToUse(j)));
            calcFRefs{j} = {grainsToUse(j),dynamicSimulatedCalcF(refPat, C.Value)};
        catch me
            disp(j)
            disp(me)
        end
%         waitbar(j/length(calcFRefs), h, "Analyzing reference patterns "+num2str(j)+"/"+num2str(length(calcFRefs)));
        send(D, 1)
    end
    saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_calcFRefs.mat");
    save(saveTo, "calcFRefs");
    close(h)
end
%%
for i = 1
    registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg";
    hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal";
    hrebsdOutputFilename = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_classic.mat");
    calcFRefsFilename = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_calcFRefs.mat");
    dicFilename = "away_screw_00"+num2str(dicSteps(i));
    reg = load(fullfile(savedir, registrationFilename)).(registrationFilename);
    hrebsdInput = load(fullfile(savedir, hrebsdFilename)).(hrebsdFilename);
    dic = load(fullfile(dicdir, dicFilename));
    calcFRefs = load(calcFRefsFilename).("calcFRefs");
%     subsetIndicies = findDataSubset(reg,dic,hrebsdInput);
    subsetIndicies = findDataSubsetFullGrains(reg,dic,hrebsdInput,minGrainSize);
%     subsetIndicies = [];
%     grains = grainsToUseFDelta{i};
%     for j = 1:length(grains)
%         subsetIndicies = [subsetIndicies; hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == grains(j))];
%     end
    hrebsdOutput = load(hrebsdOutputFilename).("hrebsdOutput");

    hrebsdOutputHybrid = hrebsdOutput;
    hrebsdOutputHybrid2 = hrebsdOutput;
    for j = 1:length(calcFRefs)
        Qps = hrebsdInput.ft.Qp2s';
        ids = hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == calcFRefs{j}{1});
        FRef = calcFRefs{j}{2}.F;
        gRef = calcFRefs{j}{2}.g;
        betaRef = calcFRefs{j}{2}.beta;
        FTest = hrebsdOutputHybrid.F(:,:,ids);
        gTest = hrebsdOutputHybrid.g(:,:,ids);
        betaTest = hrebsdOutputHybrid.beta(:,:,ids);
        [R, ~] = poldec(calcFRefs{j}{2}.F);
        hrebsdOutputHybrid2.beta(:,:,ids) = hrebsdOutput.beta(:,:,ids) - calcFRefs{j}{2}.beta;
        hrebsdOutputHybrid2.F(:,:,ids) = hrebsdOutputHybrid2.beta(:,:,ids) + eye(3);
        hrebsdOutputHybrid2.g(:,:,ids) = pagemtimes(permute(hrebsdOutputHybrid2.F(:,:,ids),[2,1,3]), calcFRefs{j}{2}.g);
        hrebsdOutputHybrid.F(:,:,ids) = pagemtimes(calcFRefs{j}{2}.F, hrebsdOutputHybrid.F(:,:,ids));
        hrebsdOutputHybrid.beta(:,:,ids) = hrebsdOutputHybrid.F(:,:,ids) - eye(3);
%         hrebsdOutputHybrid.g(:,:,ids) = pagemtimes(permute(hrebsdOutputHybrid.F(:,:,ids),[2,1,3]), calcFRefs{j}{2}.g);
    end
    post = mtexHREBSD_postProcessing(hrebsdInput, hrebsdOutputHybrid, 'DIC', dic, 'Registration', reg);
    postSubset = post.do_subset(subsetIndicies);
    post2 = mtexHREBSD_postProcessing(hrebsdInput, hrebsdOutputHybrid2, 'DIC', dic, 'Registration', reg);
    postSubset2 = post2.do_subset(subsetIndicies);
    postClassic = mtexHREBSD_postProcessing(hrebsdInput, hrebsdOutput, 'DIC', dic, 'Registration', reg);
    postClassicSubset = postClassic.do_subset(subsetIndicies);
%     saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_hybrid.mat");
%     save(saveTo, "hrebsdOutputHybrid");
%     saveTo = fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdAnalysisFilename+"_hybrid2.mat");
%     save(saveTo, "hrebsdOutputHybrid2");
end
%%
% f = fields(data);
% varPlot = "beta";
varPlot = "strain";
% varPlot = "stress";
switch varPlot
    case "beta"
        postClassicSubset.get_multiplot('beta', 'clims', [-1,1]*10E-3, 'doSmooth', 0, ...
                                       'refFrame','phosphor')
        postSubset.get_multiplot('beta', 'clims', [-1,1]*10E-3, 'doSmooth', 0, ...
                                       'refFrame','phosphor')
        postSubset2.get_multiplot('beta', 'clims', [-1,1]*10E-3, 'doSmooth', 0, ...
                                       'refFrame','phosphor')
    case "strain"
        postClassicSubset.get_multiplot('strain', 'clims', [-1,1]*4E-3, 'doSmooth', 0, ...
                                  'refFrame','phosphor')
        postSubset.get_multiplot('strain', 'clims', [-1,1]*4E-3, 'doSmooth', 0, ...
                                  'refFrame','phosphor')
        postSubset2.get_multiplot('strain', 'clims', [-1,1]*4E-3, 'doSmooth', 0, ...
                                  'refFrame','phosphor')
    case "stress"
        postClassicSubset.get_multiplot('stress', 'clims', [-1,1]*5E8, 'doSmooth', 0, ...
                                 'refFrame','dic')
        postSubset.get_multiplot('stress', 'clims', [-1,1]*5E8, 'doSmooth', 0, ...
                                 'refFrame','dic')
        postSubset2.get_multiplot('stress', 'clims', [-1,1]*5E8, 'doSmooth', 0, ...
                                 'refFrame','dic')
end

%%
clc
hrebsdOutVarName = "hrebsdOutput";%Hybrid";
hrebsdOutName = "HREBSD_gUpdate_classic";
for i = 4%:length(hrebsdSteps)
    registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg";
    hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal";
    dicFilename = "away_screw_00"+num2str(dicSteps(i));
    reg = load(fullfile(savedir, registrationFilename)).(registrationFilename);
    hrebsdInput = load(fullfile(savedir, hrebsdFilename)).(hrebsdFilename);
    dic = load(fullfile(dicdir, dicFilename));
%     subsetIndicies = findDataSubset(reg,dic,hrebsdInput);
    subsetIndicies = [];
    grains = grainsToUseFDelta{i};
    for j = 1:length(grains)
        subsetIndicies = [subsetIndicies; hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == grains(j))];
    end
    hrebsdOutput = load(fullfile(savedir, "step"+num2str(hrebsdSteps(i))+"_"+hrebsdOutName)).(hrebsdOutVarName);
    grainsToUse = unique(hrebsdInput.ebsd.grainId(subsetIndicies));
    calcFRefs = cell(length(grainsToUse), 1);
    for j = 1%1:length(calcFRefs)
        refPat = hrebsdInput.get_pattern(hrebsdInput.refIds(grainsToUse(j)));
        calcFRefs{j} = {grainsToUse(j),dynamicSimulatedCalcF(refPat, hrebsdInput)};
    end
    hrebsdOutputHybrid = hrebsdOutput;
    for j = 1%1:length(calcFRefs)
        refId = calcFRefs{j}{1};
        ids = hrebsdInput.ebsd.id(hrebsdInput.ebsd.grainId == calcFRefs{j}{1});
        hrebsdOutputHybrid.beta(:,:,ids) = hrebsdOutput.beta(:,:,ids) + calcFRefs{j}{2}.beta;%hrebsdOutputHybrid.F(:,:,ids) - eye(3);
        hrebsdOutputHybrid.F(:,:,ids) = hrebsdOutputHybrid.beta(:,:,ids) + eye(3);
        g_diff = hrebsdOutput.g(:,:,hrebsdInput.refIds(calcFRefs{j}{1})) - calcFRefs{j}{2}.g;
        hrebsdOutputHybrid.g(:,:,ids) = hrebsdOutput.g(:,:,ids) - g_diff;
        calcFRefs{j}{2}.g 
        hrebsdOutput.g(:,:,ids(4596))
        hrebsdOutputHybrid.g(:,:,ids(4596))
    end
end
%%
clc
area = [123, 141, 141, 123;
        150, 150, 165, 165];
hrebsd = load(fullfile(savedir, "test_angStep3PCCal")).("test_angStep3PCCal");
subsetIndicies = findSubset(area, hrebsd);
test_subsetFDelta = hrebsdMainFDelta(hrebsd, 'Subset', subsetIndicies);

function hrebsdMain = runHrebsdSubset(hrebsd, subsetIndicies)
    hrebsdMain = hrebsdMainClassic(hrebsd, 'Subset', subsetIndicies);
    ebsdTemp = hrebsd.ebsd(subsetIndicies);
    grainsToRun = unique(ebsdTemp.grainId);
end


function subsetIndicies = findDataSubsetFullGrains(reg, dic, hrebsd, minGrainSize)
    [x,y] = reg.transformDicCoords(dic);
    area = findArea(x,y);
    subsetIndicies2 = findSubset(area, hrebsd);
    ebsdTemp = hrebsd.ebsd(subsetIndicies2);
    ebsd = ebsdTemp(inpolygon(ebsdTemp.prop.x, ebsdTemp.prop.y, area(1,:), area(2,:)));
    grains = unique(ebsd.grainId);
    fullGrains = [];
    for i = 1:length(grains)
        ebsdGrain = hrebsd.ebsd(hrebsd.ebsd.grainId == grains(i));
        check = inpolygon(ebsdGrain.prop.x, ebsdGrain.prop.y, area(1,:), area(2,:));
        if all(check)
            if hrebsd.grains.grainSize(grains(i)) > minGrainSize
                fullGrains = [fullGrains, grains(i)];
            end
        end
    end
    subsetIndicies = [];
    for j = 1:length(fullGrains)
        subsetIndicies = [subsetIndicies; hrebsd.ebsd.id(hrebsd.ebsd.grainId == fullGrains(j))];
    end
end



function subsetIndicies = findDataSubset(reg, dic, hrebsd)
    [x,y] = reg.transformDicCoords(dic);
    area = findArea(x,y);
    subsetIndicies = findSubset(area, hrebsd);
end


function subsetIndicies = findSubset(area, hrebsd)
    ebsdArea = inpolygon(hrebsd.ebsd.prop.x, hrebsd.ebsd.prop.y, area(1,:), area(2,:));
    subsetIndicies = hrebsd.ebsd.id(ebsdArea);
end


function area = findArea(x,y)
    xMin = min(x(:));
    xMax = max(x(:));
    yMin = min(y(:));
    yMax = max(y(:));
    area = [xMin, xMax, xMax, xMin;...
            yMin, yMin, yMax, yMax];
end
