classdef mtexHREBSD_postProcessing 
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        hrebsd
        data
        plottingGrains
        F
        beta
        g
        strain
        stress
        subsetIds
        dic
        dicInterp
        registration
        Qdic
        dicData
        kappa
        alpha
        bulkGND
    end
    properties(Hidden=true)
        p
    end

    methods
        function obj = mtexHREBSD_postProcessing(mtexHREBSD,mtexHREBSD_main,varargin)
            %mtexHREBSD_postProcessing Construct an instance of this class
            %   Detailed explanation goes here
            disp("Beginning post-processing of HREBSD data...")
            obj.p = obj.create_inputParser;
            obj.hrebsd = copy(mtexHREBSD);
            obj.plottingGrains = obj.hrebsd.grains;
            obj.data = mtexHREBSD_main;
            obj.Qdic = obj.set_Qdic;
            obj.F = tensor(obj.data.F, 'rank', 2);
            obj.g = tensor(obj.data.g, 'rank', 2);
            obj.beta = obj.get_beta;
            disp("Getting strains...")
            obj.strain = obj.get_strain;
            disp("Getting stresses...")
            obj.stress = obj.get_stress;
            disp("Getting GND...")
            [obj.kappa, obj.alpha, obj.bulkGND] = obj.get_gnd;
            disp("Complete!")
            args = obj.create_additionalArgs(varargin{:});
            obj = add_args(obj, args);
        end


        function plot_bulkGND(obj, varargin)
            parse(obj.p, varargin{:});
            options = obj.p.Results;
            if options.doCbar
                [figPos, cbarPos] = obj.get_positions(options);
            else
                figPos = obj.get_positions(options);
            end
            h = figure("Units","centimeters","Position",figPos);
            if options.doSmooth
                component = obj.bulkGND;
                [ebsdGrid, newIds] = gridify(obj.hrebsd.ebsd);
                component_rs = double(ebsdGrid.isIndexed);
                component_rs(newIds) = component;
                toPlot = imgaussfilt(component_rs, 0.5);
            else
                toPlot = obj.bulkGND;
            end
            clims = obj.get_clims(toPlot(:), options);
            plot(obj.hrebsd.ebsd, toPlot, 'micronbar','off')
            set(gca, 'ColorScale', 'log')
            hold on 
            if options.doGrains && ~obj.hrebsd.isGrain
                plot(obj.plottingGrains.boundary,'lineWidth',1)
            end
            if options.RefIds
                x = obj.hrebsd.ebsd.prop.x(obj.hrebsd.refIds);
                y = obj.hrebsd.ebsd.prop.y(obj.hrebsd.refIds);
                scatter(x,y,'kx', 'LineWidth',2);
            end
            hold off
            caxis(clims)
            colormap(options.map); 
            if options.doCbar
                c = colorbar('eastoutside');
                set(get(c,'label'),'FontWeight','bold');
                set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                colormap(c, options.map)
                caxis(clims)

%                 label = obj.get_label(parameter, whichComp, options.refFrame);
%                 obj.create_colorbar(cbarPos, clims, options.map, label);
            end
            if options.doScaleBar
                [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                rectangle('Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                rectangle('Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                     'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                     'FontSize',12,'Color','w', 'FontWeight','bold')
            end
        end


        function [kappa, alpha, bulkGND] = get_gnd(obj, varargin)
            parse(obj.p, varargin{:});
            options = obj.p.Results;
            nu = options.nu;
            ebsd = obj.hrebsd.ebsd('indexed').gridify;
            kappa = ebsd.curvature;
            dS = dislocationSystem.fcc(ebsd.CS);
            dS(dS.isEdge).u = 1;
            dS(dS.isScrew).u = 1 - nu;
            dSRot = ebsd.orientations * dS;
            [rho, factor] = fitDislocationSystems(kappa,dSRot);
%             fit = dSRot.tensor .* rho;
            alpha = sum(dSRot.tensor .* rho, 2);
            bulkGND = factor*sum(abs(rho .* dSRot.u), 2);
        end


        function obj = add_args(obj, args)
            if ~isempty(args.DIC) && ~isempty(args.Registration)
                obj.dic = args.DIC;
                obj.registration = args.Registration;
                obj.dicInterp = obj.parse_DIC;
            elseif ~isempty(args.DIC) && isempty(args.Registration)
                warning("Cannont parse DIC data without registration information")
            end
        end


        function dicInterp = parse_DIC(obj)
            dicInterp = struct();
            dicComponents = {'exx', 'exy', 'eyy'};
            [X, Y] = obj.registration.transformDicCoords(obj.dic);
            for i = 1:length(dicComponents)
                comp = dicComponents{i};
                dicInterp.(comp) = scatteredInterpolant(X(:), Y(:), obj.dic.(comp)(:));
            end
        end
        

        function newObj = set_grain(obj, grainId)
            newObj = obj;
            newObj.hrebsd = copy(obj.hrebsd);
            newObj.hrebsd.isGrain = 1;
            subset = newObj.hrebsd.ebsd.grainId == grainId;
            newObj.hrebsd.ebsd = newObj.hrebsd.ebsd(subset);
            newObj.F = newObj.F(subset);
            newObj.g = newObj.g(subset);
            f = fields(newObj.strain);
            for i = 1:length(f)
                field = f{i};
                newObj.strain.(field) = newObj.strain.(field)(subset);
                newObj.stress.(field) = newObj.stress.(field)(subset);
                newObj.beta.(field) = newObj.beta.(field)(subset);
            end
        end


        function newObj = do_subset(obj, subset)
            newObj = obj;
            newObj.hrebsd = copy(obj.hrebsd);
            newObj.subsetIds = subset;
            ebsdTemp = newObj.hrebsd.ebsd(subset);
%             obj.hrebsd.ebsd = obj.hrebsd.ebsd(subset);
            xRange = [min(ebsdTemp.prop.x), max(ebsdTemp.prop.x)];
            yRange = [min(ebsdTemp.prop.y), max(ebsdTemp.prop.y)];
            xLogic = (newObj.hrebsd.ebsd.prop.x >= xRange(1)) & (newObj.hrebsd.ebsd.prop.x <= xRange(2));
            yLogic = (newObj.hrebsd.ebsd.prop.y >= yRange(1)) & (newObj.hrebsd.ebsd.prop.y <= yRange(2));
            subset2 = find(xLogic & yLogic);
            ebsdTemp2 = newObj.hrebsd.ebsd(subset2);
            [newObj.plottingGrains, ~] = calcGrains(ebsdTemp2, 'angle', 10*degree);
            newObj.hrebsd.ebsd = ebsdTemp;
            newObj.F = newObj.F(subset);
            newObj.g = newObj.g(subset);
            newObj.bulkGND = newObj.bulkGND(subset);
            newObj.alpha = newObj.alpha(subset);
            f = fields(newObj.strain);
            for i = 1:4
                field = f{i};
                newObj.strain.(field) = newObj.strain.(field)(subset);
                newObj.stress.(field) = newObj.stress.(field)(subset);
                newObj.beta.(field) = newObj.beta.(field)(subset);
            end
            [obj.hrebsd.grains, obj.hrebsd.ebsd.grainId] = calcGrains(obj.hrebsd.ebsd, 'angle',10*degree);
            if ~isempty(newObj.dic)
                newObj.dicData = struct();
                dicComponents = {'exx', 'exy', 'eyy'};
                for j = 1:length(dicComponents)
                    comp = dicComponents{j};
                    newObj.dicData.(comp) = newObj.dicInterp.(comp)(ebsdTemp.prop.x, ebsdTemp.prop.y);
                end
            end
            newObj.hrebsd.scan.Ny = length(unique(newObj.hrebsd.ebsd.prop.y));
            newObj.hrebsd.scan.Nx = length(unique(newObj.hrebsd.ebsd.prop.x));
%             [grainsLocal, ebsdLocal.grainId] = calcGrains(ebsdLocal,'angle',misorientation*degree);
        end


        function beta = get_beta(obj)
            Qps = tensor(obj.hrebsd.ft.Qp2s, 'rank', 2);
            beta = struct("phosphor", 0,...
                          "sample",0,...
                          "crystal",0,...
                          "dic", 0);
            beta.phosphor = tensor(obj.data.beta, 'rank', 2);
%             beta.sample = tensor(pagemtimes(Qps, pagemtimes(beta.phosphor.M, inv(Qps))), 'rank',2);
%             beta.sample = rotate(beta.phosphor, Qps);
%             beta.crystal = rotate(beta.sample, obj.g);

            beta.sample = Qps * beta.phosphor * Qps';
            beta.crystal = obj.g * beta.sample * obj.g';
            beta.dic = obj.Qdic * beta.sample * obj.Qdic';
        end


        function strain = get_strain(obj)
            strain = struct("phosphor", 0, ...
                            "sample", 0, ...
                            "crystal", 0, ...
                            "dic", 0);
            strain.phosphor = strainTensor(0.5*(obj.beta.phosphor + obj.beta.phosphor'));
            strain.sample = strainTensor(0.5*(obj.beta.sample + obj.beta.sample'));
            strain.crystal = strainTensor(0.5*(obj.beta.crystal + obj.beta.crystal'));
            strain.dic = strainTensor(obj.Qdic * strain.sample * obj.Qdic');
        end


        function stress = get_stress(obj)
            Qps =  tensor(obj.hrebsd.ft.Qp2s, 'rank', 2);
            stress = struct("phosphor", 0, ...
                            "sample", 0, ...
                            "crystal", 0, ...
                            "dic", 0);
            stress.crystal = obj.hrebsd.C * obj.strain.crystal*1E9; %EinsteinSum(obj.hrebsd.C, [1 2 -1 2], obj.strain.crystal, [-1 -2])
            stress.sample = stressTensor(obj.g' * stress.crystal * obj.g);
            stress.phosphor = stressTensor(Qps' * stress.sample * Qps);
            stress.dic = stressTensor(obj.Qdic * stress.sample * obj.Qdic');
        end


        function component = get_component(obj, parameter, whichComp, refFrame)
            if nargin < 3
                refFrame = 'sample';
            end
            component = obj.(parameter).(refFrame).M(whichComp(1), whichComp(2), :);
            component = component(:);
            if obj.hrebsd.isGrain
                lengthEBSD = length(obj.hrebsd.ebsd);
                component = component(1:lengthEBSD);
            end
        end


        function component_smooth = get_componentSmooth(obj, parameter, whichComp, refFrame)
            if nargin < 3
                refFrame = 'crystal';
            end
            component = obj.get_component(parameter, whichComp, refFrame);
            [ebsdGrid, newIds] = gridify(obj.hrebsd.ebsd);
            component_rs = double(ebsdGrid.isIndexed);
            component_rs(newIds) = component;
            component_filt = imgaussfilt(component_rs, 0.5);
%             component_filt = medfilt2(component_rs, [3,3]);
%             component_filt = imgaussfilt(component_filt, 0.5);
            component_smooth = component_filt(newIds);


%             component_smooth = component_smooth(newIds)
%             if ~obj.hrebsd.isGrain             
%                 component_rs = reshape(component, obj.hrebsd.scan.Ny, obj.hrebsd.scan.Nx);
%                 component_filt = imgaussfilt(component_rs, 0.75);
%                 component_smooth = component_filt(:);
%             else
%                 [ebsdGrid, newIds] = gridify(obj.hrebsd.ebsd);
%                 component_rs = double(ebsdGrid.isIndexed);
%                 component_rs(newIds) = component;
%                 component_filt = imgaussfilt(component_rs, 0.75);
%                 component_smooth = component_filt(ebsdGrid.isIndexed);
%             end
        end


        function h = plot_component(obj, parameter, whichComp, varargin)
            parse(obj.p, varargin{:});
            options = obj.p.Results;
            if options.doCbar
                [figPos, cbarPos] = obj.get_positions(options);
            else
                figPos = obj.get_positions(options);
            end
            h = figure("Units","centimeters","Position",figPos);
            if options.doSmooth
                toPlot = obj.get_componentSmooth(parameter, whichComp, options.refFrame);
            else
                toPlot = obj.get_component(parameter, whichComp, options.refFrame);
            end
            clims = obj.get_clims(toPlot(:), options);
            plot(obj.hrebsd.ebsd, toPlot, 'micronbar','off')
            hold on 
            if options.doGrains && ~obj.hrebsd.isGrain
                plot(obj.plottingGrains.boundary,'lineWidth',1)
            end
            if options.RefIds
                x = obj.hrebsd.ebsd.prop.x(obj.hrebsd.refIds);
                y = obj.hrebsd.ebsd.prop.y(obj.hrebsd.refIds);
                scatter(x,y,'kx', 'LineWidth',2);
            end
            hold off
            caxis(clims)
            colormap(options.map); 
            if options.doCbar
                c = colorbar('eastoutside');
                set(get(c,'label'),'FontWeight','bold');
                set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                colormap(c, options.map)
                caxis(clims)

%                 label = obj.get_label(parameter, whichComp, options.refFrame);
%                 obj.create_colorbar(cbarPos, clims, options.map, label);
            end
            if options.doScaleBar
                [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                rectangle('Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                rectangle('Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                     'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                     'FontSize',12,'Color','w', 'FontWeight','bold')
            end
        end


        function h = plot_componentDIC(obj, whichComp, varargin)
            parse(obj.p, varargin{:});
            options = obj.p.Results;
            dicComponents = {'exx', 'exy'; 'eyx', 'eyy'};
            comp = dicComponents{whichComp(1), whichComp(2)};
            if options.doCbar
                [figPos, cbarPos] = obj.get_positions(options);
            else
                figPos = obj.get_positions(options);
            end
            h = figure("Units","centimeters","Position",figPos);
            ax = axes(h);
%             if options.doSmooth
%                 toPlot = obj.get_componentSmooth(parameter, whichComp, options.refFrame);
%             else
%                 toPlot = obj.get_component(parameter, whichComp, options.refFrame);
%             end
            toPlot = obj.dicData.(comp);
            clims = obj.get_clims(toPlot(:), options);
            plot(obj.hrebsd.ebsd, toPlot, 'micronbar','off', 'Parent', ax)
            hold on 
            if options.doGrains && ~obj.hrebsd.isGrain
                plot(obj.plottingGrains.boundary,'lineWidth',1, 'Parent', ax)
            end
            if options.RefIds
                x = obj.hrebsd.ebsd.prop.x(obj.hrebsd.refIds);
                y = obj.hrebsd.ebsd.prop.y(obj.hrebsd.refIds);
                scatter(x,y,'kx', 'LineWidth',2);
            end
            ax.PositionConstraint = 'innerposition';
            hold off
            caxis(clims)
            colormap(options.map); 
            if options.doCbar
                c = colorbar;
                set(get(c,'label'),'FontWeight','bold');
                set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                colormap(c, options.map)
                caxis(clims)

%                 label = obj.get_label(parameter, whichComp, options.refFrame);
%                 obj.create_colorbar(cbarPos, clims, options.map, label);
            end
            if options.doScaleBar
                [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                rectangle('Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                rectangle('Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                     'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                     'FontSize',12,'Color','w', 'FontWeight','bold')
            end
        end


        function multiplotDIC(obj, varargin)
            parse(obj.p, varargin{:})
            options = obj.p.Results;
            figPos = obj.get_positions(options);
            label = "\epsilon [\mum/\mum]";
            fig = figure("Units","centimeters","Position",figPos);
            tiles = tiledlayout(fig, 2, 2, "TileSpacing","compact", "Padding","compact");
            dicComponents = {'exx', 'exy', 'cbar', 'eyy'};
            labels = {'xx', 'xy', 'cbar','yy'};
            for i = 1:length(dicComponents)
                comp = dicComponents{i};
                ax = nexttile(tiles);
                if i == 3
                    axis off
                    c = colorbar(ax, 'north');
                    set(get(c,'label'),'string', label,'FontWeight','bold');
                    set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                    colormap(c, options.map)
                    caxis(clims)
                else
                    hold on
                    toPlot = obj.dicData.(comp);
                    clims = obj.get_clims(toPlot(:), options);
                    plot(obj.hrebsd.ebsd, toPlot, 'parent', ax)
                    if options.doGrains && ~obj.hrebsd.isGrain
                        plot(obj.plottingGrains.boundary,'lineWidth',1,'parent', ax)            
                    end
                    pos = get(ax, "Position");
                    labelPos = [pos(1) + pos(3)/2, pos(2),0,0];
                    annotation('textbox', labelPos, 'string', labels{i}, ...
                        'HorizontalAlignment','center','FontName','Times New Roman', ...
                        'FontSize',12)
                    axis off
                    hold off
                    xlabel(labels(i));
                    caxis(clims);
                    colormap(options.map);
                end      
                if i == 4 
                    [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                    rectangle('Parent',ax,'Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                    rectangle('Parent',ax,'Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                    text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                         'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                         'FontSize',12,'Color','w', 'FontWeight','bold')
                end
            end
        end



        function get_multiplot(obj, parameter, varargin)
            parse(obj.p, varargin{:})
            options = obj.p.Results;
            switch options.Multiplot
                case 'all'
                    obj.multiplotFull(parameter, options)
                case '2d'
                    obj.multiplot2D(parameter, options)
            end
        end


        function multiplot2D(obj, parameter, options)
            figPos = obj.get_positions(options);
            label = obj.get_parameterLabel(parameter);
            fig = figure("Units","centimeters","Position",figPos);
            tiles = tiledlayout(fig, 2, 2, "TileSpacing","compact", "Padding","compact");
%             labels = {'xx', 'xy', 'cbar','yy'};
            if options.refFrame == "dic"
                labels = {'xx', 'xy'; 'cbar','yy'};
            else
                labels = {'11', '12'; 'cbar','22'};
            end
            for i = 1:2
                for j = 1:2
                    ax = nexttile(tiles);
                    axPos = get(ax, "Position");
                    if i == 2 && j == 1
                        axis off
                        c = colorbar(ax, 'north');
                        if parameter == "stress"
                            set(get(c,'label'),'string', "\sigma [GPa]",'FontWeight','bold');
                            set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                            colormap(c, options.map)
                            caxis(clims/1E9)
                        else
                            set(get(c,'label'),'string', label+"_{"+options.refFrame+"}",'FontWeight','bold');
                            set(c,'FontName', 'Times New Roman', 'FontSize', 12)
                            colormap(c, options.map)
                            caxis(clims)
                        end
                    else
                        hold on 
                        if options.doSmooth
                            toPlot = obj.get_componentSmooth(parameter, [i,j], options.refFrame);
                        else
                            toPlot = obj.get_component(parameter, [i,j], options.refFrame);
                        end
%                         labels{i,j}
%                         text(ax,0.025, 0.6, labels{i,j})
                        clims = obj.get_clims(toPlot(:), options);
                        plot(obj.hrebsd.ebsd, toPlot, 'parent', ax)
                        if options.doGrains && ~obj.hrebsd.isGrain
                            plot(obj.plottingGrains.boundary,'lineWidth',1,'parent', ax)            
                        end
                        if options.RefIds
                            refIds = obj.get_refIds;
                            x = obj.hrebsd.ebsd.prop.x(refIds);
                            y = obj.hrebsd.ebsd.prop.y(refIds);
                            scatter(x,y,'kx', 'LineWidth',2,'parent', ax);
                        end
                        pos = get(ax, "Position");
                        labelPos = [pos(1) + pos(3)/2, pos(2),0,0];
                        annotation('textbox', labelPos, 'string', labels{i,j}, ...
                            'HorizontalAlignment','center','FontName','Times New Roman', ...
                            'FontSize',12)
                        axis off
                        hold off
                        caxis(clims);
                        colormap(options.map);
                    end
                    if i == 2 && j == 2
                        [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                        rectangle('Parent',ax,'Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                        rectangle('Parent',ax,'Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                        text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                             'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                             'FontSize',12,'Color','w', 'FontWeight','bold')
                    end
                end
            end
%             c = colorbar;
%             c.Layout.Tile = 'east';
%             set(get(c,'label'),'string', label+"_{"+options.refFrame+"}");
%             set(c,'FontName', 'Times New Roman', 'FontSize', 12)
%             colormap(c, options.map)
%             caxis(clims)
        end


        function multiplotFull(obj, parameter, options)
            figPos = obj.get_positions(options);
            label = obj.get_parameterLabel(parameter);
            fig = figure("Units","centimeters","Position",figPos);
            tiles = tiledlayout(fig, 3, 3, "TileSpacing","compact", "Padding","compact");
            for i = 1:3
                for j = 1:3
                    ax = nexttile(tiles);
                    hold on 
                    if options.doSmooth
                        toPlot = obj.get_componentSmooth(parameter, [i,j], options.refFrame);
                    else
                        toPlot = obj.get_component(parameter, [i,j], options.refFrame);
                    end
                    clims = obj.get_clims(toPlot(:), options);
                    plot(obj.hrebsd.ebsd, toPlot, 'parent', ax)
                    if options.doGrains && ~obj.hrebsd.isGrain
                        plot(obj.plottingGrains.boundary,'lineWidth',1,'parent', ax)            
                    end
                    if options.RefIds
                        grainsToUse = unique(obj.hrebsd.ebsd.grainId);
                        refIds = obj.hrebsd.refIds(grainsToUse);
%                         refIds = obj.get_refIds(grainsToUse);
                        x = obj.hrebsd.ebsd.prop.x(obj.hrebsd.ebsd.id == refIds);
                        y = obj.hrebsd.ebsd.prop.y(refIds);
                        scatter(x,y,'kx', 'LineWidth',2,'parent', ax);
                    end
                    axis off
                    hold off
                    caxis(clims);
                    colormap(options.map);
                end
                if i == 3 && j == 3
                    [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                    rectangle('Parent',ax,'Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                    rectangle('Parent',ax,'Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                    text(textPos(1), textPos(2), num2str(sbSize)+"\mum", ...
                         'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                         'FontSize',12,'Color','w', 'FontWeight','bold')
                end
            end
            c = colorbar;
            c.Layout.Tile = 'east';
            set(get(c,'label'),'string', label+"_{"+options.refFrame+"}");
            set(c,'FontName', 'Times New Roman', 'FontSize', 12)
            colormap(c, options.map)
            caxis(clims)
        end


        function refIds = get_refIds(obj)
            if obj.hrebsd.isGrain
                refIds = obj.hrebsd.ebsd.id == obj.hrebsd.refIds;
            else
                refIds = obj.hrebsd.refIds;
            end
        end


        function strainVec = tensor2vector(obj, tensor, H, numpoints)
            if tensor.size(1) == 1
                strainVec = [tensor.M(1,1), tensor.M(2,2), tensor.M(3,3),...
                             tensor.M(2,3), tensor.M(1,3), tensor.M(1,2)]';
            else
                strainVec = zeros(6, tensor.size(1));
                for i = 1:tensor.size(1)
                    strainVec(:, i) = obj.tensor2vector(tensor(i));
                    if nargin > 2 && mod(i,30) == 0
                        percentComplete = i/(numpoints/2)*100;
                        waitbar(i/numpoints, H, "Converting strain tensors to vectors... " + num2str(percentComplete,3) + "%")
                    end
                end
            end
        end


        function tensor = vector2tensor(obj, vector, H, numpoints)
            if size(vector,2) == 1
                tensor = [vector(1), vector(6), vector(5);...
                          vector(6), vector(2), vector(4);...
                          vector(5), vector(4), vector(3)];
            else
                tensor = zeros(3,3,size(vector,2));
                for i = 1:size(vector,2)
                    tensor(:,:,i) = obj.vector2tensor(vector(:,i));
                    if nargin > 2 && mod(i,30) == 0
                        percentComplete = i/(numpoints/2)*100;
                        waitbar((i+numpoints/2)/numpoints, H, "Converting stress vectors to tensors..." + num2str(percentComplete,3) + "%")
                    end
                end
            end
        end


        function [posOuter, posInner, posText, scalebarSize] = get_scalebarPosition(obj)
            xrange = [min(obj.hrebsd.ebsd.x), max(obj.hrebsd.ebsd.x)];
            yrange = [min(obj.hrebsd.ebsd.y), max(obj.hrebsd.ebsd.y)];
            hfw = max(obj.hrebsd.ebsd.x) - min(obj.hrebsd.ebsd.x);
            vfw = max(obj.hrebsd.ebsd.y) - min(obj.hrebsd.ebsd.y);
            scalebarSize = obj.get_scalebarSize(hfw);
            sbHeight = 0.125*vfw;
            pad = 0.025*hfw;
            posOuter = [xrange(1)+1,yrange(2)-sbHeight-1,scalebarSize+2*pad, sbHeight];
            posInner = [xrange(1)+1+pad, yrange(2)-sbHeight+sbHeight*3/5, scalebarSize, sbHeight/4];
            posText = [(xrange(1)+1+pad)+(scalebarSize)/2, yrange(2)-sbHeight*0.8];
        end
    end


    methods(Static)
        function scalebarSize = get_scalebarSize(hfw)
            if hfw > 10
                scalebarSize = ceil(0.25*hfw/10)*10;
            end
        end

%         function dicInterp = parse_DIC(dic, reg)
%             [X, Y] = reg.transformDicCoords(dic);
%             keys = ['exx', 'eyy', 'exy'];
%             dicInterp = struct;
%             for i = 1:length(keys)
%                 s.(keys(i)) = struct;
%                 s.(keys(i)).f = scatteredInterpolant(X(:), Y(:), dic.(keys(i))(:));
%             end
%         end


        function Qdic = set_Qdic
            c = cos(-pi/2);
            s = sin(-pi/2);
            Qdic = tensor([c, -s, 0;
                           s,  c, 0;
                           0,  0, 1], ...
                           'rank', 2);
        end


        function label = get_parameterLabel(parameter)
            switch parameter
                case "strain"
                    label = "\epsilon";
                case "beta"
                    label = "\beta";
                case "stress"
                    label = "\sigma";
            end
        end


        function label = get_label(parameter, whichComp, refFrame)
            i = whichComp(1);
            j = whichComp(2);
            switch parameter
                case "strain"
                    parameter_label = "\epsilon";
                case "beta"
                    parameter_label = "\beta";
                case "stress"
                    parameter_label = "\sigma";
            end
            label = string(refFrame) + parameter_label+"_{"+num2str(i)+num2str(j)+"}";
        end

        function create_colorbar(cbarPos, clims, map, label)
            h = colorbar; 
            set(h,'FontName', 'Times New Roman', 'FontSize', 14)
            set(h, 'units', 'centimeters', 'Position', cbarPos)
            set(get(h,'label'),'string', label);
            caxis(clims); 
            colormap(map); 
        end

        function p = create_inputParser
            p = inputParser;
            addParameter(p, 'map', jet(256));
            addParameter(p, 'clims', 2);
            addParameter(p, 'figSize', [18.3,18.3]);
            addParameter(p, 'doCbar', 1);
            addParameter(p, 'doScaleBar', 1);
            addParameter(p, 'doGrains', 1);
            addParameter(p, 'doSmooth', 1);
            addParameter(p, 'RefIds', 0);
            addParameter(p, 'refFrame', 'crystal');
            addParameter(p, 'Multiplot', 'all')
            addParameter(p, 'nu', 0.3)
        end


        function args = create_additionalArgs(varargin)
            parser = inputParser;
            addParameter(parser, 'DIC', []);
            addParameter(parser, 'Registration', []);
            parse(parser, varargin{:});
            args = parser.Results;
        end


        function [figPos, varargout] = get_positions(options)
            if options.doCbar
                cbarPos = [options.figSize(1)-0.2*options.figSize(1),...
                           0.1*options.figSize(2),...
                           options.figSize(1)*0.02,...
                           options.figSize(1) - 0.35*options.figSize(2)];
                figPos = [5,5,options.figSize(1),...
                    options.figSize(2)-0.15*options.figSize(2)];
                if nargout == 2
                    varargout = {cbarPos};
                end
            else
                figPos = [5,5,options.figSize(1),...
                    options.figSize(2)-0.15*options.figSize(2)];
            end
        end


        function clims = get_clims(A, options)
            % Default is to use clims with +/- 2std
            if length(options.clims) == 1
                meanA = mean(A);
                stdA = std(A);
                clims = [meanA - options.clims*stdA, meanA + options.clims*stdA];
            else
                clims = options.clims;
            end
        end
    end
end