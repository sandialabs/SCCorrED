classdef mtexHREBSD_postProcessing < dynamicprops
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        hrebsd
        data
        plottingGrains
        omega
        subsetIds
        dic
        dicInterp
        registration
        Qdic
        dicData
    end
    properties(Hidden=true)
        p
    end

    methods
        function obj = mtexHREBSD_postProcessing(mtexHREBSD,mtexHREBSDmain,varargin)
            %mtexHREBSD_postProcessing Construct an instance of this class
            %   Detailed explanation goes here
            disp("Beginning post-processing of HREBSD data...")
            obj.p = obj.create_inputParser;
            obj.hrebsd = copy(mtexHREBSD);
            obj.plottingGrains = obj.hrebsd.grains;
            obj.data = mtexHREBSDmain;
            obj.Qdic = obj.set_Qdic;
            if isprop(mtexHREBSDmain, 'beta')
                checkAddProp(obj, 'F');
                checkAddProp(obj, 'g');
                checkAddProp(obj, 'beta');
                checkAddProp(obj, 'strain');
                checkAddProp(obj, 'stress');
                obj.F = tensor(obj.data.F, 'rank', 2);
                obj.g = tensor(obj.data.g, 'rank', 2);
                obj.beta = obj.get_beta;
                disp("Getting strains...")
                obj.strain = obj.get_strain;
                disp("Getting stresses...")
                obj.stress = obj.get_stress;
            end
            if isprop(mtexHREBSDmain, 'dbetadx')
                Qps = tensor(obj.hrebsd.ft.Qp2s, 'rank', 2);
%                 beta.sample = Qps * beta.phosphor * Qps';
                checkAddProp(obj, 'dbetadx');
                checkAddProp(obj, 'dbetady');
                ax = tensor(mtexHREBSDmain.dbetadx/(mtexHREBSD.scan.xStep*1e-6), 'rank', 2);
                bx = Qps * ax * Qps';
                ay = tensor(mtexHREBSDmain.dbetady/(mtexHREBSD.scan.yStep*1e-6), 'rank', 2);
                by = Qps * ay * Qps';
                obj.dbetadx = bx.M;
                obj.dbetady = by.M;
                checkAddProp(obj, 'alpha');
                checkAddProp(obj, 'alphaSmooth');
                checkAddProp(obj, 'bulkGND')
                checkAddProp(obj, 'bulkGNDSmooth')
                obj.getAlpha;
            end

%             disp("Getting GND...")
%             [obj.kappa, obj.alpha, obj.bulkGND] = obj.get_gnd;
            disp("Complete!")
            args = obj.create_additionalArgs(varargin{:});
            obj = add_args(obj, args);
            if isprop(mtexHREBSDmain, 'Subset')
                obj = obj.do_subset(mtexHREBSDmain.Subset);
            end
        end


        function obj = getAlpha(obj)
            [ebsdGrid, ~] = gridify(obj.hrebsd.ebsd);
            dxgrid = reshape( ...
                obj.dbetadx, [3,3,size(ebsdGrid,2), size(ebsdGrid,1)] ...
                );
            dygrid = reshape( ...
                obj.dbetady, [3,3,size(ebsdGrid,2), size(ebsdGrid,1)] ...
                );
            alphaGrid = obj.partialcurl(dxgrid, dygrid);
%             alphaGrid = obj.partialcurl(obj.dbetadx, obj.dbetady);
            obj.alpha = reshape(alphaGrid, [6,length(obj.hrebsd.ebsd)]);
%             obj.alpha = obj.partialcurl(obj.dbetadx, obj.dbetady);
            dxsmooth = zeros(size(dxgrid));
            dysmooth = zeros(size(dygrid));
            for i = 1:9
                j = ceil(i/3);
                k = mod(i-1,3) + 1;
                dxsmooth(j,k,:,:) = imgaussfilt(squeeze(dxgrid(j,k,:,:)), 0.5);
                dysmooth(j,k,:,:) = imgaussfilt(squeeze(dygrid(j,k,:,:)), 0.5);
            end
            alphaSmoothGrid = obj.partialcurl(dxsmooth, dysmooth);
            obj.alphaSmooth = reshape( ...
                alphaSmoothGrid, [6,length(obj.hrebsd.ebsd)] ...
                );
            alphasum = sum(abs(alphaGrid),1);
%             alphasum = sum(abs(obj.alpha),1);
            gndgrid = alphasum/(287*1e-12);
            obj.bulkGND = reshape(gndgrid, size(obj.hrebsd.ebsd));
            alphasumsmooth = sum(abs(alphaSmoothGrid),1);
            gndgridsmooth = alphasumsmooth/(287*1e-12); 
            obj.bulkGNDSmooth = reshape(gndgridsmooth, size(obj.hrebsd.ebsd));
        end


        function obj = alphaPost(obj)
            X = obj.hrebsd.ebsd.prop.x;
            Y = obj.hrebsd.ebsd.prop.y;
            nx = obj.hrebsd.scan.Nx;
            ny = obj.hrebsd.scan.Ny;
            alphalocal = zeros([6,size(obj.data.beta,3)]);
            obj.bulkGND = zeros([1,size(obj.data.beta,3)]);
            xstep = obj.hrebsd.scan.xStep;
            ystep = obj.hrebsd.scan.yStep;
            x = 1;
            y = 1;
            for i = 1:size(alphalocal,2)
                if x == nx
                    x = 1;
                    y = y + 1;
                    continue
                elseif y == ny
                    break
                else
                    id = sub2ind([nx,ny], x, y);
                    idx = sub2ind([nx,ny], x+1, y);
                    idy = sub2ind([nx,ny], x, y+1);
    
                    b = obj.beta.sample.M(:,:,id);
                    bx = obj.beta.sample.M(:,:,idx);
                    by = obj.beta.sample.M(:,:,idy);
                    dbx = (bx - b)/(xstep*1e-6);
                    dby = (by - b)/(ystep*1e-6);
                    alphalocal(:,i) = obj.partialcurl(dbx, dby);
                    x = x + 1;
                end
            end
            obj.alpha = alphalocal;
            alphasum = sum(abs(obj.alpha),1);
            obj.bulkGND = alphasum/(287*1e-12);
        end 


        function plot_bulkGND(obj, varargin)
            parse(obj.p, varargin{:});
            options = obj.p.Results;
            if options.doCbar
                [figPos, cbarPos] = obj.get_positions(options);
            else
                figPos = obj.get_positions(options);
            end
            h = figure;
            if options.doSmooth
                toPlot = obj.bulkGNDSmooth;
%                 [ebsdGrid, newIds] = gridify(obj.hrebsd.ebsd);
%                 component_rs = double(ebsdGrid.isIndexed);
%                 component_rs(newIds) = component;
%                 toPlot = imgaussfilt(component_rs, 0.5);
            else
                toPlot = obj.bulkGND;
            end
            % remove values less than or equal to zero
            toPlot(toPlot <= 0) = 1e-3;
            clims = obj.get_clims(toPlot(:), options);
            clims = [13,16];
            plot( ...
                obj.hrebsd.ebsd, log10(toPlot), ...
                'micronbar','off')
%             set(gca, 'ColorScale', 'log')
            hold on 
            if options.doGrains && ~obj.hrebsd.isGrain
                plot( ...
                    obj.plottingGrains.boundary, ...
                    'lineWidth',1, ...
                    'lineColor', options.lineColor);
            end
            if options.RefIds
                x = obj.hrebsd.ebsd.prop.x(obj.ebsd.id == obj.hrebsd.refIds);
                y = obj.hrebsd.ebsd.prop.y(obj.ebsd.id ==obj.hrebsd.refIds);
                scatter(x,y,'kx', 'LineWidth',2);
            end
            hold off
            caxis(clims)
            colormap(options.map); 
            if options.doCbar
                c = mtexColorbar('eastoutside');
                set(get(c,'label'),'FontWeight','bold', 'String', "log_{10}(bulk GND [m^{-2}])");
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
                text(textPos(1), textPos(2), num2str(sbSize)+" \mum", ...
                     'HorizontalAlignment', 'center', 'FontName','Times New Roman',...
                     'FontSize',10,'Color','w', 'FontWeight','bold')
            end

            % set figure size
            set(h,'Units','centimeters','Position',figPos);
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
            if isprop(obj, 'beta')
                newObj.hrebsd.ebsd = ebsdTemp;
                newObj.F = newObj.F(subset);
                newObj.g = newObj.g(subset);
    %             newObj.bulkGND = newObj.bulkGND(subset);
    %             newObj.alpha = newObj.alpha(subset);
                f = fields(newObj.strain);
                for i = 1:4
                    field = f{i};
                    newObj.strain.(field) = newObj.strain.(field)(subset);
                    newObj.beta.(field) = newObj.beta.(field)(subset);
                    if ~isempty(newObj.stress)
                        newObj.stress.(field) = newObj.stress.(field)(subset);
                    end
                end
            end
            if isprop(obj, 'alpha')
                newObj.alpha = newObj.alpha(:, subset);
                newObj.alphaSmooth = newObj.alphaSmooth(:, subset);
                newObj.dbetadx = newObj.dbetadx(:,:, subset);
                newObj.dbetady = newObj.dbetady(:,:, subset);
                newObj.bulkGND = newObj.bulkGND(subset);
                newObj.bulkGNDSmooth = newObj.bulkGNDSmooth(subset);
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
        end


        function beta = get_beta(obj)
            Qps = tensor(obj.hrebsd.ft.Qp2s, 'rank', 2);
            beta = struct("phosphor", 0,...
                          "sample",0,...
                          "crystal",0,...
                          "dic", 0);
            beta.phosphor = tensor(obj.data.beta, 'rank', 2);
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
            stress = struct("phosphor", [], ...
                            "sample", [], ...
                            "crystal", [], ...
                            "dic", []);
            stress.crystal = stressTensor;
            stress.phosphor = stressTensor;
            stress.sample = stressTensor;
            stress.dic = stressTensor;

            % reassign values with correct stiffness tensor

            for i = 1:size(obj.hrebsd.C, 2)
                if ~isempty(obj.hrebsd.C{i})
                    C = obj.hrebsd.C{i};
                    phase = C.CS.mineral;
                    ids = obj.hrebsd.ebsd(phase).id;
                    stress.crystal(ids) = ...
                       C * obj.strain.crystal(ids) * 1E9; % Pa
                    stress.sample(ids) = stressTensor( ...
                        obj.g(ids)' * stress.crystal(ids) * obj.g(ids) ...
                        );
                    stress.phosphor(ids) = stressTensor( ...
                        Qps' * stress.sample(ids) * Qps ...
                        );
                    stress.dic(ids) = stressTensor( ...
                        obj.Qdic * stress.sample(ids) * obj.Qdic' ...
                        );
                end
            end

%             stress.crystal = obj.hrebsd.C{2} * obj.strain.crystal * 1E9;
%             stress.sample = stressTensor(obj.g' * stress.crystal * obj.g);
%             stress.phosphor = stressTensor(Qps' * stress.sample * Qps);
%             stress.dic = stressTensor(obj.Qdic * stress.sample * obj.Qdic');

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


        function refIds = get_refIds(obj)
            if obj.hrebsd.isGrain
                refIds = obj.hrebsd.ebsd.id == obj.hrebsd.refIds;
            else
                refIds = obj.hrebsd.refIds;
            end
        end


%         function strainVec = tensor2vector(obj, tensor, H, numpoints)
%             if tensor.size(1) == 1
%                 strainVec = [tensor.M(1,1), tensor.M(2,2), tensor.M(3,3),...
%                              tensor.M(2,3), tensor.M(1,3), tensor.M(1,2)]';
%             else
%                 strainVec = zeros(6, tensor.size(1));
%                 for i = 1:tensor.size(1)
%                     strainVec(:, i) = obj.tensor2vector(tensor(i));
%                     if nargin > 2 && mod(i,30) == 0
%                         percentComplete = i/(numpoints/2)*100;
%                         waitbar(i/numpoints, H, "Converting strain tensors to vectors... " + num2str(percentComplete,3) + "%")
%                     end
%                 end
%             end
%         end


%         function tensor = vector2tensor(obj, vector, H, numpoints)
%             if size(vector,2) == 1
%                 tensor = [vector(1), vector(6), vector(5);...
%                           vector(6), vector(2), vector(4);...
%                           vector(5), vector(4), vector(3)];
%             else
%                 tensor = zeros(3,3,size(vector,2));
%                 for i = 1:size(vector,2)
%                     tensor(:,:,i) = obj.vector2tensor(vector(:,i));
%                     if nargin > 2 && mod(i,30) == 0
%                         percentComplete = i/(numpoints/2)*100;
%                         waitbar((i+numpoints/2)/numpoints, H, "Converting stress vectors to tensors..." + num2str(percentComplete,3) + "%")
%                     end
%                 end
%             end
%         end


        function [posOuter, posInner, posText, scalebarSize] = get_scalebarPosition(obj)
            xrange = [min(obj.hrebsd.ebsd.x), max(obj.hrebsd.ebsd.x)];
            yrange = [min(obj.hrebsd.ebsd.y), max(obj.hrebsd.ebsd.y)];
            hfw = max(obj.hrebsd.ebsd.x) - min(obj.hrebsd.ebsd.x);
            vfw = max(obj.hrebsd.ebsd.y) - min(obj.hrebsd.ebsd.y);
            % Automatic scalebar sizing
            sbFrac = 0.25;
            expTerm = 10^floor(log10(sbFrac*hfw));
            possibleSizes = [1,2.5,5];
            [~,i] = min(abs((sbFrac*hfw/expTerm./possibleSizes - 1)));
            scalebarSize = possibleSizes(i) * expTerm;
            sbHeight = 0.125*vfw;
            pad = 0.02*hfw;
            posOuter = [xrange(1)+pad,yrange(2)-sbHeight-pad,scalebarSize+2*pad, sbHeight];
            posInner = [xrange(1)+2*pad, yrange(2)-sbHeight*2/5-pad, scalebarSize, sbHeight/4];
            posText = [(xrange(1)+2*pad)+(scalebarSize)/2, yrange(2)-sbHeight*0.8-pad];
        end


        function checkAddProp(obj, propname)
            if ~isprop(obj, propname)
                addprop(obj, propname);
            end
        end
    end


    methods(Static)
        function alpha = partialcurl(betaderiv1,betaderiv2)
            alpha = zeros(6,size(betaderiv1,3),size(betaderiv1,4));
            
            alpha(1,:,:) = betaderiv2(1,1,:,:) - betaderiv1(1,2,:,:); %13
            alpha(2,:,:) = betaderiv2(2,1,:,:) - betaderiv1(2,2,:,:); %23
            alpha(3,:,:) = betaderiv2(3,1,:,:) - betaderiv1(3,2,:,:); %33
            alpha(4,:,:) = betaderiv1(1,3,:,:); %12
            alpha(5,:,:) = -betaderiv2(2,3,:,:); %21
            alpha(6,:,:) = -betaderiv2(1,3,:,:) - betaderiv1(2,3,:,:); %11 - 22
        end


        function id = getId(X, Y, x, y)
            xlogic = X == x;
            ylogic = Y == y;
            id = find(xlogic.*ylogic);
        end

        function Qdic = set_Qdic
%             c = cos(-pi/2);
%             s = sin(-pi/2);
%             Qdic = tensor([c, -s, 0;
%                            s,  c, 0;
%                            0,  0, 1], ...
%                            'rank', 2);
            c = cos(-pi);
            s = sin(-pi);
            Qdic = tensor([1, 0,  0;
                           0,  c,-s;
                           0,  s, c], ...
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
                case "omega"
                    label = "\omega";
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
                case "omega"
                    parameter_label = "\omega";
            end
            label = string(refFrame) + parameter_label+"_{"+num2str(i)+num2str(j)+"}";
        end

        function create_colorbar(cbarPos, clims, map, label)
            h = colorbar; 
            set(h,'FontName', 'Times New Roman', 'FontSize', 10)
            set(h, 'units', 'centimeters', 'Position', cbarPos)
            set(get(h,'label'),'string', label);
            clim(clims); 
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
            addParameter(p, 'lineColor', 'k')
            fig1 = figure('visible','off'); % Fix for preventing empty figure
            addParameter(p, 'Parent', gca)
        end


        function args = create_additionalArgs(varargin)
            parser = inputParser;
            addParameter(parser, 'DIC', []);
            addParameter(parser, 'Registration', []);
            addParameter(parser, 'SplitDislocation', 0);
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


        obj = multiplotStrainBeta(obj, varargin)
    end
end