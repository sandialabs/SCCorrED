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
            % Apply mask if supplied
            if isfield(options,"noMaskIds")
                indicies = ~ismember( ...
                    obj.hrebsd.ebsd.id, options.noMaskIds ...
                    );
                toPlot(indicies) = nan;
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
                y = obj.hrebsd.ebsd.prop.y(obj.hrebsd.ebsd.id == refIds);
                scatter(x,y,'kx', 'LineWidth',2,'parent', ax);
            end
            axis off
            hold off
            caxis(clims);
            colormap(options.map);
        end
        if i == 3 && j == 3
            if options.doScaleBar
                [pos1, pos2, textPos, sbSize] = obj.get_scalebarPosition;
                rectangle('Parent',ax,'Position',pos1, 'EdgeColor','k', 'FaceColor','k')
                rectangle('Parent',ax,'Position',pos2, 'EdgeColor','w', 'FaceColor','w')
                text(textPos(1), textPos(2), ...
                     strcat(num2str(sbSize)," ",char(956),'m'), ...
                     'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                     'FontSize',10,'Color','w', 'FontWeight','bold')
            end
        end
    end
    switch parameter
        case 'strain'
            units = "\epsilon";
        case 'beta'
            units = "\epsilon";
        case 'stress'
            units = 'Pa';
    end
    c = colorbar;
    c.Layout.Tile = 'south';
    set(get(c,'label'),'string', label+"_{"+options.refFrame+"} ["+units+"]");
    set(c,'FontName', 'Times New Roman', 'FontSize', 12)
    c.FontSize = 12;
    c.Ruler.Exponent = 0;
%     cbar.tick_params(labelsize=10) 
    colormap(c, options.map)
    caxis(clims)
end