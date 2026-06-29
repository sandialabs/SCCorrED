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