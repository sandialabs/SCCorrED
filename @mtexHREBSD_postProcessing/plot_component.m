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