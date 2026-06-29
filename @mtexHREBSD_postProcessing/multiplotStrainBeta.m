function multiplotStrainBeta(obj, varargin)
    options = parseKwargs(varargin{:});
    options = set_clims(options);
    figPos = obj.get_positions(options);
    fig = figure("Units","centimeters","Position",figPos);
    tiles = tiledlayout( ...
        fig, 10, 3, "TileSpacing","compact", "Padding","compact" ...
        );
%     climStrain = get_clim()
    for i = 1:3
        for j = 1:3
            ax = nexttile(tiles, [3,1]);
            toPlot = get_vals(obj, [i,j], options);
            if options.noMaskIds
                indicies = ~ismember( ...
                    obj.hrebsd.ebsd.id, options.noMaskIds ...
                    );
                toPlot(indicies) = nan;
            end
            plot(obj.hrebsd.ebsd, toPlot, 'parent', ax, 'micronbar','off')
            if options.RefIds
                hold on
                grainsToUse = unique(obj.hrebsd.ebsd.grainId);
                refIds = obj.hrebsd.refIds(grainsToUse);
%                         refIds = obj.get_refIds(grainsToUse);
                x = obj.hrebsd.ebsd.prop.x(obj.hrebsd.ebsd.id == refIds);
                y = obj.hrebsd.ebsd.prop.y(obj.hrebsd.ebsd.id == refIds);
                scatter(x,y,'kx', 'LineWidth',2,'parent', ax);
            end
            if options.doGrains && ~obj.hrebsd.isGrain
                plot(obj.plottingGrains.boundary,'lineWidth',1,'parent', ax)            
            end
            ax.Position;
            caxis(get_clim(toPlot, options.("clims"+isBetaComp([i,j]))))
            colormap(ax, options.map)
            t = text( ...
                ax, 0, 0,get_label([i,j]), ...
                'VerticalAlignment','top', ...
                'HorizontalAlignment','left', ...
                'Color', options.labelColor,...
                'FontSize',options.FontSize, ...
                'FontWeight', options.FontWeight, ...
                'FontName',options.FontName,...
                'BackgroundColor','k', ...
                'Margin',1 ...
                );
            set(t, 'Units','normalized')
            set(t, 'Position', [0,1])
            if i == 3
                if j == 1
                    pos31 = ax.Position;
                elseif j == 3
                    pos33 = ax.Position;
                end
            end
        end
    end
    ax = nexttile(tiles);
    pos = ax.Position;
    set(ax, 'Visible', 'off')
    c1 = colorbar(ax,'Location', 'south', ...
        'Position', [pos31(1),pos(2)+pos(4)-0.015,pos31(3),0.015] ...
        );
    set(get(c1,'label'),'string',"\omega");
    set(c1,'FontName', options.FontName, 'FontSize', options.FontSize)
    caxis(options.climsBeta)
    cbar_handle = findobj(gcf,'tag','Colorbar');
    set(cbar_handle, 'XAxisLocation','bottom')
    colormap(ax, options.map)

    ax = nexttile(tiles); % skip middle tile

    set(ax, 'Visible', 'off')
    ax = nexttile(tiles);
    pos = ax.Position;
    set(ax, 'Visible', 'off')
    c2 = colorbar(ax,'Location', 'south', ...
        'Position', [pos33(1),pos(2)+pos(4)-0.015,pos33(3),0.015] ...
        );
    set(get(c2,'label'),'string',"\epsilon");
    set(c2,'FontName', options.FontName, 'FontSize', options.FontSize)
    caxis(options.climsStrain)
    cbar_handle = findobj(gcf,'tag','Colorbar');
    set(cbar_handle(1), 'XAxisLocation','bottom')
    colormap(ax, options.map)
end


function options = set_clims(options)
    if ~isempty(options.clims)
        options.climsStrain = options.clims;
        options.climsBeta = options.clims;
    end
end


function clims = get_clim(A, opt)
    % Default is to use clims with +/- 2std
    if length(opt) == 1
        meanA = mean(A);
        stdA = std(A);
        clims = [meanA - opt*stdA, meanA + opt*stdA];
    else
        clims = opt;
    end
end


function field = isBetaComp(comp)
    check = ismember(comp, [2,1;3,1;3,2], 'rows');
    if check
        field = 'Beta';
    else
        field = 'Strain';
    end
end


function label = get_label(comp)
    if ismember(comp, [2,1;3,1;3,2], 'rows')
        label= "\omega_{"+num2str(comp(1))+num2str(comp(2)) + "}";
    else
        label = "\epsilon_{"+num2str(comp(1))+num2str(comp(2)) + "}";
    end
end


function toPlot = get_vals(obj, comp, options)
    if ismember(comp, [2,1;3,1;3,2], 'rows')
        parameter = 'beta';
        if options.doSmooth
            A = obj.get_componentSmooth( ...
                parameter, comp, options.refFrame ...
                );
            B = obj.get_componentSmooth( ...
                'strain', comp, options.refFrame ...
                );
        else
            A = obj.get_component(parameter, comp, options.refFrame);
            B = obj.get_component('strain', comp, options.refFrame);
        end
        toPlot = A - B;
    else
        parameter = 'strain';
        if options.doSmooth
            toPlot = obj.get_componentSmooth( ...
                parameter, comp, options.refFrame ...
                );
        else
            toPlot = obj.get_component(parameter, comp, options.refFrame);
        end
    end
end


function options = parseKwargs(varargin)
    p = inputParser;
    addParameter(p, 'map', jet(256));
    addParameter(p, 'noMaskIds', []);
    addParameter(p, 'clims', []);
    addParameter(p, 'climsStrain', 2);
    addParameter(p, 'climsBeta', 2);
    addParameter(p, 'figSize', [18.3,18.3]);
    addParameter(p, 'doCbar', 1);
    addParameter(p, 'doScaleBar', 1);
    addParameter(p, 'doGrains', 1);
    addParameter(p, 'doSmooth', 1);
    addParameter(p, 'RefIds', 0);
    addParameter(p, 'refFrame', 'crystal');
    addParameter(p, 'labelColor', 'white');
    addParameter(p, 'FontSize', 12);
    addParameter(p, 'FontWeight', 'normal');
    addParameter(p, 'FontName', "Times New Roman");
    parse(p, varargin{:});
    options = p.Results;
end