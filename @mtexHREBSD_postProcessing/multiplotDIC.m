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