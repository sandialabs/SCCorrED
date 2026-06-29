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
                % Apply mask if supplied
                if isfield(options,"noMaskIds")
                    indicies = ~ismember( ...
                        obj.hrebsd.ebsd.id, options.noMaskIds ...
                        );
                    toPlot(indicies) = nan;
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
    end