clc; clear; %close all;
hrebsdSteps = [0,1,2,3];
dicSteps = [2,3,4,5];

figDir = 'C:\Users\wggilli\Documents\Papers\SEM_ DIC_and_HREBSD_data_set\Figures';
hrebsdOutVarName = "hrebsdOutput";
hrebsdOutName = "HREBSD_gUpdate_classic";
loaddir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting';
savedir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting\PDR_PCCal';
%%
% dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Away_screw";
% step0Grains = [438, 741, 748, 467, 456, 516 ,700, 589, 515, 733, 775, 724, 621, 613, 565, 548, 604, 729, 562, 651, 599, 639, 491];
% fnamePrefix = "farScrew_";

dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Near_screw";
step0Grains = [760,613,685,646,648,929,889,777,841,587,589,585,881]; %593,597
fnamePrefix = "nearScrew_";
%%
data = struct();
for i = 1:length(hrebsdSteps)
%     "near_screwStep0PCCal_updatedRefs2"
    if i == 1
        hrebsdFilename = "near_screwStep"+num2str(hrebsdSteps(i))+"PCCal_updatedRefs2";
%         hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal_fullTractionFree_updatedRefs";
        hrebsdInput = load(fullfile(loaddir, hrebsdFilename)).("test_angStep0PCCal");
    else
%         hrebsdFilename = "test_angStep"+num2str(hrebsdSteps(i))+"PCCal";
        hrebsdFilename = "near_screwStep"+num2str(hrebsdSteps(i))+"PCCal";
        hrebsdInput = load(fullfile(loaddir, hrebsdFilename)).("test_angStep"+num2str(hrebsdSteps(i)+"PCCal"));
    end
%     hrebsdInput = load(fullfile(loaddir, hrebsdFilename)).(hrebsdFilename);
    key = "Step"+num2str(hrebsdSteps(i));
    data.(key) = struct('inputs', hrebsdInput);
end

data.Step0.grains = step0Grains;
%%
for i = 1:length(data.Step0.grains)
    grainInt = data.Step0.grains(i);
    x1 = data.Step0.inputs.grains(grainInt).centroid;
%     figure;
%     plot(data.Step0.inputs.grains(grainInt))  
    for j = 2:length(hrebsdSteps)
        key = "Step"+num2str(hrebsdSteps(j));
        x2 = data.(key).inputs.grains.centroid;
        dist = ((x2(:,1) - x1(1)).^2 + (x2(:,2) - x1(2)).^2).^0.5;
%         [~, grainStep] = min(dist)
        diffGrainSize = abs(data.Step0.inputs.grains(grainInt).area - data.(key).inputs.grains.area);
        distN = (dist - min(dist))/(max(dist) - min(dist));
        diffGSN = (diffGrainSize - min(diffGrainSize))/(max(diffGrainSize) - min(diffGrainSize));
        [~, grainStep] = min(distN+diffGSN);
        data.(key).grains(i) = grainStep;
%         figure
%         plot(data.(key).inputs.grains(grainStep))
    end
end

%%
clc
f = fields(data);
totalPoints = 0;
for i = 1:length(f)
    key = f{i};
    data.(key).subset = [];
    for j = 1:length(data.(key).grains)
        data.(key).subset = [data.(key).subset; data.(key).inputs.ebsd.id(data.(key).inputs.ebsd.grainId == data.(key).grains(j))];
    end
    totalPoints = totalPoints + length(data.(key).subset);
end
% refPat = data.Step0.inputs.get_pattern(data.Step0.inputs.refIds(data.Step0.grains(1)));
refPats = cell(1, length(data.Step0.grains));
for i = 1:length(refPats)
    refPats{i} = data.Step0.inputs.get_pattern(data.Step0.inputs.refIds(data.Step0.grains(i)));
end
%%
clc
R = parallel.pool.Constant(refPats);
% inputs = {data.Step0.inputs, data.Step1.inputs, data.Step2.inputs, data.Step3.inputs}
C = parallel.pool.Constant(data);
outputs = cell(3,1);
h = waitbar(0, "Starting PDR analysis");
nUpdateWaitbar(totalPoints, h);
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateWaitbar)
tic
parfor i = 1:length(f)
    key = f{i};
    a = C.Value.(key).inputs;
    b = mtexHREBSD_main(a, 1);
    for k = 1:length(C.Value.(key).grains)
        refPat = R.Value{k};
        grains = C.Value.(key).grains;
        ids = a.ebsd.id(a.ebsd.grainId == grains(k));
        for j = 1:length(ids)
            testPat = a.get_pattern(ids(j));
            try
                calcF = classicCalcF(refPat, testPat, a);
        %         disp(calcF.F)
                b.beta(:,:,ids(j)) = calcF.beta;
        %         disp(calcF.beta)
                b.F(:,:,ids(j)) = calcF.F;
                b.g(:,:,ids(j)) = calcF.g;
                b.fitMetrics.SSE(ids(j)) = calcF.fit.metrics.SSE;
                b.fitMetrics.R2(ids(j)) = calcF.fit.metrics.R2;
            catch me 
                disp(me)
            end
            send(D, 1)
        end
    end
    outputs{i} = b;
end

for i = 1:length(f)
    key = f{i};
    data.(key).output = outputs{i};
end
saveTo = fullfile(savedir, fnamePrefix + "PDR_Data_updatedRefs2");
save(saveTo, "data");
%%
clc

dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Away_screw";
step0Grains = [438, 741, 748, 467, 456, 516 ,700, 589, 515, 733, 775, 724, 621, 613, 565, 548, 604, 729, 562, 651, 599, 639, 491];
fnamePrefix = "farScrew_";
data = load(fullfile(savedir, "PDR_Data_fullTractionFree_updatedRefs")).("data");

% dicdir = "\\snl\home\wggilli\Documents\2023\Microstructure_clones\Dec_29_2022\images\Near_screw";
% step0Grains = [760,613,685,646,648,929,889,777,841,587,589,585,881]; %593,597
% fnamePrefix = "nearScrew_";
% data = load(fullfile(savedir, fnamePrefix + "PDR_Data_fullTractionFree_updatedRefs_nearScrew")).("data");
% 
% data.Step0.grains = step0Grains;
%%
for i = 1:length(hrebsdSteps)
    registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg";
    dicFilename = "away_screw_00"+num2str(dicSteps(i));
%     registrationFilename = "step"+num2str(hrebsdSteps(i))+"Reg_nearScrew";
%     dicFilename = "near_screw_00"+num2str(dicSteps(i));
    reg = load(fullfile(loaddir, registrationFilename)).(registrationFilename);
    dic = load(fullfile(dicdir, dicFilename));
    key = "Step"+num2str(hrebsdSteps(i));
%     data.(key).output = outputs{i};
    data.(key).reg = reg;
    data.(key).dic = dic;
    data.(key).subset = [];
    for j = 1:length(data.(key).grains)
        data.(key).subset = [data.(key).subset; data.(key).inputs.ebsd.id(data.(key).inputs.ebsd.grainId == data.(key).grains(j))];
    end
    data.(key).post = mtexHREBSD_postProcessing(data.(key).inputs,data.(key).output, 'DIC', dic, 'Registration', reg);
    data.(key).postSubset = data.(key).post.do_subset(data.(key).subset);
end

%%
% figDir = 'C:\Users\wggilli\Desktop\ToTimM&M';
clc
f = fields(data);
% varPlot = "beta";
% varPlot = "strain";
varPlot = "stress";
for i = 1:length(f)
    switch varPlot
        case "beta"
            data.(f{i}).postSubset.get_multiplot('beta', 'clims', [-1,1]*5E-3, 'doSmooth', 0, ...
                                           'refFrame','phosphor')
%             data.(f{i}).postSubset.get_multiplot('beta', 'clims', [-1,1]*2E-3, 'doSmooth', 0, ...
%                                            'refFrame','crystal')
        case "strain"
            data.(f{i}).postSubset.get_multiplot('strain', 'clims', [-1,1]*1E-2, 'doSmooth', 0, ...
                                      'refFrame','dic', 'Multiplot', '2d')
            data.(f{i}).postSubset.multiplotDIC('clims', [-1,1]*1E-2)
        case "stress"
            data.(f{i}).postSubset.get_multiplot('stress', 'clims', [-1,1]*1.5E9, 'doSmooth', 1, ...
                                     'refFrame','dic', 'Multiplot', '2d')
%             saveTo = fullfile(figDir, fnamePrefix+"PDRstress_"+f{i}+".png");
%             exportgraphics(gcf,saveTo,'BackgroundColor','white')
%             close(gcf);
%             data.(f{i}).postSubset.multiplotDIC('clims', [-1,1]*3E-2)
%             saveTo = fullfile(figDir, fnamePrefix+"strain_"+f{i}+".png");
%             exportgraphics(gcf,saveTo,'BackgroundColor','white')
%             close(gcf);
    end
end

%% Figures for first paper
close all
f = fields(data);
for i = 1:length(f)
    if i == length(f)
        doScaleBar = 1;
    else
        doScaleBar = 0;
    end
%     data.(f{i}).postSubset.plot_component('stress', [2,2], 'clims', [-1,1]*1.5E9, ...
%         'doSmooth', 1, 'refFrame', 'dic', 'doScaleBar', doScaleBar, ...
%         'figSize', [5,5], 'doCbar', 0)
%     saveTo = fullfile(figDir, fnamePrefix+"PDRstress_22_"+f{i}+".png");
%     exportgraphics(gcf,saveTo,'BackgroundColor','white')
%     close(gcf);

    data.(f{i}).postSubset.plot_componentDIC([2,2], 'clims', [-1,1]*3E-2,...
        'doScaleBar', doScaleBar, 'figSize', [5,5], 'doCbar', 0);
%     saveTo = fullfile(figDir, fnamePrefix+"DICstrain_22_"+f{i}+".png");
%     exportgraphics(gcf,saveTo,'BackgroundColor','white')
%     close(gcf);
end
%%
close all; clc
f = fields(data);
fig = figure("Units","centimeters","Position",[5,5,18.9,8]);
tiledlayout(fig,8,17,"TileSpacing","compact", "Padding","compact")
for i = 1:length(f)
    ax = nexttile([4,4]);
    toPlot = data.(f{i}).postSubset.dicData.('eyy');
    plot(data.(f{i}).postSubset.hrebsd.ebsd, toPlot, 'Parent', ax)
    colormap(ax, 'jet')
    caxis(ax, [-1,1]*3E-2)
    hold on 
    plot(data.(f{i}).postSubset.plottingGrains.boundary,'lineWidth',1, 'Parent', ax)
    pos = get(ax, "Position");
    labelPos = [pos(1) + pos(3)/2, 0.5,0,0];
    title(f{i},'FontName','Times New Roman','FontSize',12,'Fontweight','normal',...
        'Units', 'normalized', 'Position', [0.5, -0.25, 0])
%     annotation('textbox', labelPos, 'string', f{i}, ...
%         'HorizontalAlignment','center','VerticalAlignment','middle', ...
%         'FontName','Times New Roman','FontSize',12)
    
    if i == 1
        pos = get(ax, "Position");
        ht = annotation('textarrow', [pos(1)-0.02,pos(1)-0.02], [pos(2)+pos(4)/2,pos(2)+pos(4)/2],...
            'string', 'Strain','HeadStyle', 'none', 'LineStyle', 'none', ...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontName','Times New Roman','FontSize',12,'TextRotation',90);
    end
    axis off
    hold off
    if i == length(f)
        ax = nexttile([4,1]);
        c = colorbar(ax,'east');
        set(get(c,'label'),'string', "\epsilon [\mum/\mum]");
%         set(c.XLabel,{'String','Rotation','Position'},{"\epsilon",0,[0.5 0.045]})
        set(c, 'YAxisLocation','right')
        set(c,'FontName', 'Times New Roman', 'FontSize', 12)
        colormap(c, 'jet')
        caxis([-1,1]*3E-2)
        c.Ruler.Exponent = -2;
        axis off
    end
end

for i = 1:length(f)
    ax = nexttile([4,4]);
    toPlot = data.(f{i}).postSubset.get_componentSmooth('stress', [2,2], 'dic');
    plot(data.(f{i}).postSubset.hrebsd.ebsd, toPlot, 'Parent', ax)
    colormap(ax, 'jet')
    caxis(ax,[-1,1]*1.5E9)
    hold on 
    plot(data.(f{i}).postSubset.plottingGrains.boundary,'lineWidth',1, 'Parent', ax)
%     title(f{i},'FontName','Times New Roman','FontSize',12,'Fontweight','normal')
    if i == 1
        pos = get(ax, "Position");
        ht = annotation('textarrow', [pos(1)-0.02,pos(1)-0.02], [pos(2)+pos(4)/2,pos(2)+pos(4)/2],...
            'string', 'Stress','HeadStyle', 'none', 'LineStyle', 'none', ...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontName','Times New Roman','FontSize',12,'TextRotation',90);
    end
    axis off
    hold off
    if i == length(f)
        ax = nexttile([4,1]);
        c = colorbar('east');
        set(get(c,'label'),'string', "\sigma [GPa]");
%         set(c.XLabel,{'String','Rotation','Position'},{"\sigma",0,[0.5 0.5]})
        set(c, 'YAxisLocation','right')
        set(c,'FontName', 'Times New Roman', 'FontSize', 12)
        colormap(c, 'jet')
        caxis([-1,1]*1.5)
        axis off
    end
end
saveTo = fullfile(figDir, fnamePrefix+"stressStrainComp.png");
exportgraphics(gcf,saveTo,'BackgroundColor','white')
close(gcf);

%%
close all
fig = figure("Units","centimeters","Position",[5,5,1,5])
fig.Units = 'normalized';
axis off
c = colorbar('west');
set(gca,'Visible', false)
c.Position = [0.5, 0.15, 0.15, 0.74]
% c.Label.String = "\epsilon"
set(get(c,'label'),'String', "\epsilon",'FontWeight','bold');
set(c,'FontName', 'Times New Roman', 'FontSize', 12)
colormap(c, 'jet')
caxis([-1,1]*3E-2)
%%
% clc
% dicComps = {'exx', 'exy', 'eyy'};
% hrebsdComps = {[1,1], [1,2], [2,2]};
% s = struct('exx', struct('stress', [], 'strain', []),...
%            'exy', struct('stress', [], 'strain', []),...
%            'eyy', struct('stress', [], 'strain', []));
% for i = 1:length(f)
%     dicComp = dicComps;
%     key = f{i};
%     for j = 1:length(data.(key).grains)
%         for k = 1:length(dicComps)
%             dicComp = dicComps{k};
%             l = hrebsdComps{k};
%             ids = find(data.(key).postSubset.hrebsd.ebsd.grainId == data.(key).grains(j));
% %         inds = data.(key).postSubset.hrebsd.ebsd.id == ids
%             strain = mean(data.(key).postSubset.dicData.(dicComp)(ids));
%             stress = mean(data.(key).postSubset.stress.dic.M(l(1), l(2), ids));
% %             s.(dicComp).strain{i} = strain;
% %             s.(dicComp).stress{i} = stress;
%             s.(dicComp).strain = [s.(dicComp).strain; strain];
%             s.(dicComp).stress = [s.(dicComp).stress; stress(:)];
% %             plot(strain, stress);
% %         disp(data.(key).grains(j))
%         end
%     end
% end
% %%
% close all
% dicComps = {'exx', 'exy', 'dummy', 'eyy'};
% ylims = [[-8E8,2E8]; [-1,1]*6E8; [-1,1]; [-0.5,1]*1.5E9]/1E6;
% fig = figure;
% tiledlayout(fig, 2, 2, "TileSpacing","compact", "Padding","compact");
% for i = 1:length(dicComps)
%     ax = nexttile;
%     if i == 3
%         axis off
%     else
%         c = dicComps{i};
% %         if i == 1
% %             scatter(abs(s.(c).strain), abs(s.(c).stress)/1E6, 250, 'b.')
% %         else
%         scatter(s.(c).strain, s.(c).stress/1E6, 200, 'b.')
% %         end
%         grid on
%         ylim(ylims(i,:))
%         xlabel('\epsilon')
%         ylabel('\sigma [MPa]')
%     end
% end
% %%
% f = fields(data);
% for i = 1:length(f)
%     datai = data.(f{i});
%     figure; plot(datai.inputs.grains);
% end
% 
% %%
% ebsd = data.Step3.postSubset.hrebsd.ebsd;
% [ebsdGrid, newId] = gridify(ebsd)
% %%
% grains = [146, 146, 151, 143];
% refPat = data.Step0.inputs.get_pattern(data.Step0.inputs.refIds(grains(1)))
% outputs = {};
% C = parallel.pool.Constant(data);
% parfor i = 1:length(f)
%     a = C.Value.(f{i}).inputs;
%     b = mtexHREBSD_main(a, 1)
%     beta = zeros(3,3,totalPoints);
%     g = zeros(3,3,totalPoints);
%     ids = a.ebsd.id(a.ebsd.grainId == grains(i));
%     for j = 1:length(ids)
%         testPat = a.get_pattern(ids(j));
%         calcF = classicCalcF(refPat, testPat, a);
%         b.beta(:,:,ids(j)) = calcF.beta;
%         b.F(:,:,ids(j)) = calcF.F;
%         b.g(:,:,ids(j)) = calcF.g;
%         b.fitMetrics.SSE(ids(j)) = calcF.fit.metrics.SSE;
%         b.fitMetrics.R2(ids(j)) = calcF.fit.metrics.R2;
%     end
%     outputs{i} = b;
% end