classdef mtexHREBSD_calcF2 < dynamicprops
    %mtexHREBSD_calcF is the super class for calcF methods
    %   Detailed explanation goes here

    %TODO: Move mtexHREBSD_pattern SimData property to mtexHREBSD_calcF2, 
    % require simulation analysis type to initialize

    properties
        ft 
        roi
        assumption {mustBeMember(...
            assumption,{'free-surface', 'traction-free','trace=0'}...
            )} = 'free-surface'
        AnalysisType (1,1) {...
                    mustBeInteger,...
                    mustBeMember(AnalysisType, [0,1,2,3])...
                    } = 0
        RobustFit (1,1) {...
                    mustBeInteger,...
                    mustBeMember(RobustFit, [0,1])...
                    } = 0
        roiCenters_p
        Tolerance
        IterationLimit = 20;
%         numStds = 1.5
        MADThreshold = 0.025
    end


    properties(Transient, Hidden = true)
        cctiles
        useq
    end


    methods
        function obj = mtexHREBSD_calcF2(input, kwargs)
            %mtexHREBSD_calcF Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                input mtexHREBSD
                kwargs.RobustFit = 0
                kwargs.AnalysisType = 0
                kwargs.Assumption = []
                kwargs.Tolerance (1,1) double = input.analysis.stepTolerance;
%                 kwargs.numStds = []
%                 kwargs.stdThreshold = []
                kwargs.MADThreshold = [];
                kwargs.Simulators {...
                    mustBeA(...
                        kwargs.Simulators,...
                        ["EBSPSim","double"]...
                    )} = []
            end
%             disp(kwargs)
%             obj = obj.parseKWargs(varargin{:});
            obj.RobustFit = kwargs.RobustFit;
            obj.AnalysisType = kwargs.AnalysisType;
            obj.Tolerance = kwargs.Tolerance;
%             if kwargs.numStds 
%                 obj.numStds = kwargs.numStds;
%             end
            if kwargs.MADThreshold
                obj.MADThreshold = kwargs.MADThreshold;
            end
            if ~isempty(kwargs.Simulators)
                obj.checkAddProp('Simulators');
                obj.Simulators = kwargs.Simulators;
            end
%             obj.residual = mtexHREBSD_fit("SSE", "R2");
            obj.ft = input.ft;
            obj.roi = input.roi;
            obj.roiCenters_p = obj.ft.imageVec2phosphorFrame( ...
                obj.roi.centers ...
            );
            if isempty(kwargs.Assumption)
                obj.assumption = input.analysis.assumptions;
            else
                obj.assumption = kwargs.Assumption;
            end
        end


        function obj = getBeta(obj, testPattern, refPattern, options)
            arguments
                obj
                testPattern mtexHREBSD_pattern
                refPattern {mustBeA(refPattern, {'mtexHREBSD_pattern', 'double'})} = []
                options.FGuess (3,3) double = eye(3)
                options.plotCC (1,1) {...
                    mustBeMember(options.plotCC, [1,0])...
                    } = 0
                options.Verbose (1,1) {...
                    mustBeMember(options.Verbose, [1,0])...
                    } = 0 
                options.ConvergencePlots (1,1) {...
                    mustBeMember(options.ConvergencePlots, [1,0])...
                    } = 0 
                options.Tolerance (1,1) double = 0
                options.NormXCorr (1,1) {...
                    mustBeMember(options.NormXCorr, [1,0])...
                    } = 0 
                options.MakeGIF = []
                options.GIFfps (1,1) double = 3;
                options.SaveIterations (1,1) {...
                    mustBeMember(options.SaveIterations, [1,0])...
                    } = 0 
            end
            if ~options.Tolerance
                options.Tolerance = obj.Tolerance;
            end
            if options.SaveIterations
                obj.checkAddProp('ImageIterations')
                obj.checkAddProp('qIteration')
                obj.qIteration = zeros( ...
                    3, obj.roi.numRois, obj.IterationLimit +1 ...
                    );
                obj.ImageIterations = zeros( ...
                    obj.roi.pixelSize, ...
                    obj.roi.pixelSize, ...
                    obj.IterationLimit +1 ...
                    );
                obj.ImageIterations(:,:,1) = ...
                    uint8(rescale(testPattern.image, 0, 255));
            end
            obj.checkAddProp('beta');
            obj.checkAddProp('q');
            obj.checkAddProp('qFit');
            obj.checkAddProp('qFitU');
            obj.checkAddProp('qFitR');
            obj.checkAddProp('g');
            obj.checkAddProp('R2');
            obj.checkAddProp('SSE');
            if obj.AnalysisType == 0
                if isempty(refPattern)
                    error("Classic analysis requires a reference pattern!")
                else
%                     disp(testPattern)
%                     disp(refPattern)
                    obj = obj.betaClassic(testPattern, refPattern, options);
                end
            elseif obj.AnalysisType == 1
                obj = obj.betaDynSim(testPattern, options);
            elseif obj.AnalysisType == 2
                obj = obj.betaFDelta(testPattern, options);
            elseif obj.AnalysisType == 3
                obj.checkAddProp('CC');
                obj = obj.betaXASGO(refPattern, testPattern, options);
            end
        end


        function g = getg(obj, g)
            F = obj.beta + eye(3);
            [r,~] = poldec(obj.ft.Qp2s*F*obj.ft.Qp2s');
            g = g*r';
        end


        function checkAddProp(obj, propname)
            if ~isprop(obj, propname)
                addprop(obj, propname);
            end
        end


        obj = betaClassic(obj, testPattern, refPattern, options);
        obj = betaDynSim(obj, testPat, options);
        obj = betaFDelta(obj, testPattern, options);
        

        function out = calc_deviatoricVR(obj)
            [V,R] = polarDecompositionVR(obj.beta + eye(3,3));
            I = eye(3,3);
            out = (I + V - trace(V)*I./3.0)*R;
        end


        function obj = calc_betaDetector( ...
                obj, refPattern, rs, qs ...
                )
            r1 = rs(1,obj.useq)';
            r2 = rs(2,obj.useq)';
            r3 = rs(3,obj.useq)';
            q1 = qs(1,obj.useq)';
            q2 = qs(2,obj.useq)';
            [A4, b4] = obj.construct_Ab(refPattern,rs,qs);
            if obj.RobustFit
                [X,~] = robustfit(A4,b4,'bisquare',4.685,'off');
            else
                X = A4\b4;
            end
            obj.beta = reshape(X, [3 3])';
            F = obj.beta + eye(3);
            [R,U] = poldec(F);
            strain = U - eye(3);
            rot = R - eye(3);
            strainT = strain'; 
            RT = rot';
            bFitU = A4 * strainT(:);
            bFitR = A4 * RT(:);
            bFit = A4*X;
            q1Fit = bFit(1:length(q1))./r3;
            q1FitU = bFitU(1:length(q1))./r3;
            q1FitR = bFitR(1:length(q1))./r3;
            fitInd = length(q1)*2;
            q2Fit = bFit(length(q1)+1:fitInd)./r3;
            q2FitU = bFitU(length(q1)+1:fitInd)./r3;
            q2FitR = bFitR(length(q1)+1:fitInd)./r3;
            qFit = [q1Fit, q2Fit]';
            qFitU = [q1FitU, q2FitU]';
            qFitR = [q1FitR, q2FitR]';
            qFit(3,:) = 0;
            qFitU(3,:) = 0;
            qFitR(3,:) = 0;
            obj.qFit = qFit;
            obj.qFitU = qFitU;
            obj.qFitR = qFitR;
            obj.R2 = obj.get_R2(obj.q(1:2,obj.useq), qFit(1:2,:));
            residual = ((obj.q(:,obj.useq) - qFit).^2);
            obj.SSE = obj.get_SSE(residual);
        end


        function plot(obj, pattern, varargin)
            p = inputParser;
            addParameter(p, 'multiplier',1)
            addParameter(p, 'showFit',0)
            addParameter(p, 'showRotStretch',0)
            addParameter(p, 'showMeas_Rot', 0)
            addParameter(p, 'numRois', obj.roi.numRois)
            addParameter(p, 'showRois', 0)
            addParameter(p, 'roiSize', 0)
            addParameter(p, 'Parent', gcf)
            parse(p, varargin{:});
            options = p.Results;
            m = options.multiplier;
            pattern.show_pattern;
            scaledImage = uint8(rescale(pattern.image,0,255));
            imshow(scaledImage);
            colormap gray
            hold on
%             colors = ['r', 'b'];
%             for i = 1:2
%             if i == 1
%             q_iv = obj.ft.phosphorFrame2imageVec(obj.q(:,~obj.useq));
%             centersPlot = obj.roi.centers(:,~obj.useq);
%             else
                q_iv = obj.ft.phosphorFrame2imageVec(obj.q(:,logical(obj.useq)));
                centersPlot = obj.roi.centers(:,logical(obj.useq));
%             end
            numSkip = round(obj.roi.numRois/options.numRois);
            centers = obj.roi.centers(:, 2:numSkip:obj.roi.numRois);
            qPlot = q_iv;%(:, 2:numSkip:obj.roi.numRois);
            q = quiver( ...
                centersPlot(1,:), centersPlot(2,:), ...
                m*qPlot(1,:), ...
                m*qPlot(2,:), ...
                'b', 'LineWidth',1, ...
                'DisplayName', 'Measured' ...
                );
            q.AutoScale = 'off';
%             end
            if options.showFit
                qfit_iv = obj.ft.phosphorFrame2imageVec(obj.qFit);
                qFitU_iv = obj.ft.phosphorFrame2imageVec(obj.qFitU);
                qfitPlot = qfit_iv(:, 2:numSkip:obj.roi.numRois);
                qfitUPlot = qFitU_iv(:, 2:numSkip:obj.roi.numRois);
                qfit = quiver( ...
                    centers(1,:), centers(2,:), ...
                    m*qfitPlot(1,:), ...
                    m*qfitPlot(2,:), ...
                    'g','LineWidth',1, ...
                    'DisplayName', 'Fit' ...
                    );
                qfit.AutoScale = 'off';
                % qfitU = quiver( ...
                %     centers(1,:), centers(2,:), ...
                %     m*qfitUPlot(1,:), ...
                %     m*qfitUPlot(2,:), ...
                %     'r','LineWidth', 1, ...
                %     'Displayname', 'FitU' ...
                %     );
                % qfitU.AutoScale = 'off';
                legend
            end
%             qPlot = qPlot - qfitPlot;
%             q = quiver( ...
%                 centers(1,:), centers(2,:), ...
%                 m*qPlot(1,:), ...
%                 m*qPlot(2,:), ...
%                 'b', 'LineWidth',1, ...
%                 'DisplayName', 'Measured' ...
%                 );
%             q.AutoScale = 'off';
            if options.showRotStretch
                qfitR = obj.ft.phosphorFrame2imageVec(obj.qFitR);
                qfitU = obj.ft.phosphorFrame2imageVec(obj.qFitU);
                centersR = centers + m.*qfitU(:, 2:numSkip:obj.roi.numRois);
                
                qfitPlot = qfitR(:, 2:numSkip:obj.roi.numRois);
                qfitR = quiver( ...
                    centersR(1,:), centersR(2,:), ...
                    m*qfitPlot(1,:), ...
                    m*qfitPlot(2,:), ...
                    'y','LineWidth',1, ...
                    'DisplayName', 'Rotation' ...
                    );
                qfitR.AutoScale = 'off';
                
                qfitPlot = qfitU(:, 2:numSkip:obj.roi.numRois);
                qfitU = quiver( ...
                    centers(1,:), centers(2,:), ...
                    m*qfitPlot(1,:), ...
                    m*qfitPlot(2,:), ...
                    'r','LineWidth',1, ...
                    'DisplayName', 'Stretch' ...
                    );
                qfitU.AutoScale = 'off';
                legend
            end
            if options.showMeas_Rot
                qfitR = obj.ft.phosphorFrame2imageVec(obj.qFitR);
                qfitU = obj.ft.phosphorFrame2imageVec(obj.qFitU);
                % q_iv = obj.ft.phosphorFrame2imageVec(obj.q(:,logical(obj.useq)));
                q_meas_rot = q_iv - qfitR;                               
                
                qfitPlot = q_meas_rot(:, 2:numSkip:obj.roi.numRois);
                q_meas_rot = quiver( ...
                    centers(1,:), centers(2,:), ...
                    m*qfitPlot(1,:), ...
                    m*qfitPlot(2,:), ...
                    'c','LineWidth',1, ...
                    'DisplayName', 'Meas-Rot' ...
                    );
                q_meas_rot.AutoScale = 'off';

                qfitPlot = qfitU(:, 2:numSkip:obj.roi.numRois);
                qfitU = quiver( ...
                    centers(1,:), centers(2,:), ...
                    m*qfitPlot(1,:), ...
                    m*qfitPlot(2,:), ...
                    'r','LineWidth',1, ...
                    'DisplayName', 'Stretch' ...
                    );
                qfitU.AutoScale = 'off';
                legend
            end
            if options.roiSize
                for i = 1:length(centers)
                    x = [centers(1,i)-options.roiSize/2; centers(1,i)+options.roiSize/2];
                    y = [centers(2,i)-options.roiSize/2; centers(2,i)+options.roiSize/2];
                    plot([x(1), x(2), x(2), x(1), x(1)], [y(1), y(1), y(2), y(2), y(1)], 'Color', 'white')
                end
            end
            hold off
            axis equal
        end


        function plotCC(obj, ccMap, testroi, refroi, q)
            if size(q,1) == 2
                q_iv = q;
            else
                q_iv = obj.ft.phosphorFrame2imageVec(q);    
            end
%             q_iv = q;
            figure;
            obj.cctiles = tiledlayout(1,3);
            ax = nexttile(obj.cctiles);
            hold(ax,'on')
            scaledImage = uint8(rescale(testroi,0,255));
            imshow(scaledImage);
            axis(ax, 'equal')
            title("q_x = " + num2str(q_iv(1)))
            ax = nexttile(obj.cctiles);
            scaledImage = uint8(rescale(refroi,0,255));
            imshow(scaledImage);
            hold(ax,'on')
            axis(ax, 'equal')
            title("q_y = " + num2str(q_iv(2)))
            ax = nexttile(obj.cctiles);
            imagesc(ax, ccMap)
            colorbar
            axis(ax, 'image')
            colormap(ax, jet(256))
%             axis off
            [~,I] = max(ccMap, [], 'all');
            [y,x] = ind2sub(size(ccMap), I);
            title("Peak at (" +num2str(x) + "," + num2str(y) + ")")
            xlabel(ax,"dx = " +num2str(size(ccMap,1)/2 - x))
            ylabel(ax,"dy = " +num2str(size(ccMap,2)/2 - y))
        end


        function movie(obj, n, fps, kwargs)
            arguments
                obj mtexHREBSD_calcF2
                n double = 1 %TODO add vector functionality
                fps double = 5
                kwargs.mult (1,1) double = 1
                kwargs.SaveToGIF = []
                kwargs.PlayMovie = 1
            end
            if ~isprop(obj, 'ImageIterations')
                error(['No iteration data stored in "' inputname(1) '". ' ...
                    'Use Name/Value pair" ''SaveIterations'', 1 " ' ...
                    'when calling mtexHREBSD_calcF2.getBeta'])
            end
            if length(n) > 1
                arange = n(2:end);
            else
                arange = 1:size(obj.ImageIterations,3);
            end
            frames = struct('cdata',[],'colormap',[]); % Why does it make me do this matlab?
            frames(length(arange)) = ...
                struct('cdata',[],'colormap',[]);
            centersPlot = obj.roi.centers;
            figure;
            for i = 1:length(arange)
%                 frames(i).cdata = ...
%                     repmat(uint8(obj.ImageIterations(:,:,i)), [1,1,3]);
%                 figure;
                j = arange(i);
                imshow(uint8(obj.ImageIterations(:,:,j)), 'Parent',gca)
                qPlot = obj.ft.phosphorFrame2imageVec( ...
                    obj.qIteration(:,:,j) ...
                    );
                hold on
                q = quiver( ...
                    centersPlot(1,:), centersPlot(2,:), ...
                    kwargs.mult*qPlot(1,:), ...
                    kwargs.mult*qPlot(2,:), ...
                    'r', 'LineWidth',1, ...
                    'DisplayName', 'Measured' ...
                    );
                q.AutoScale = 'off';
                hold off
                xlabel(j)
                frames(i) = getframe(gcf);
                pause(1/fps)
%                 close gcf
            end
            if ~isempty(kwargs.SaveToGIF)
                obj.makeGIF(frames, kwargs.SaveToGIF, fps);
            end
            if kwargs.PlayMovie == 1
                movie(gcf, frames, n(1), fps);
            end
        end
    end

    methods(Static)

%         function goodIndicies = getGoodIndicies(qs)
%             [~, outlierInd] = rmoutliers(qs');
%             goodIndicies = ~outlierInd';
%         end


        % Updated by TJB, 2024-11-22
        function R2 = get_R2(a, b)
            cov_mat = cov(a,b); % Covariance matrix
            % Pearson's correlation coefficient
            pcc = cov_mat(1,2) / sqrt(cov_mat(1,1) * cov_mat(2,2));
            R2 = pcc^2;
        end


        function SSE = get_SSE(residual)
            SSE = sum(residual(:)).^0.5;
        end


        function makeGIF(frames, fname, fps)
            filename = fname;
            for i = 1:length(frames)
                [A, map] = rgb2ind(frames(i).cdata, 256);
                if i == 1
                    imwrite(A, map, filename,"gif", ...
                        "LoopCount",Inf, ...
                        "DelayTime",1/fps);
                else
                    imwrite(A, map,filename,"gif", ...
                        "WriteMode","append", ...
                        "DelayTime",1/fps);
                end
            end
            [path, file, ext] = fileparts(fname);
            filename = fullfile( ...
                string(path), string(file)+"Converged"+string(ext) ...
                );
            for i = [1,length(frames)]
                [A, map] = rgb2ind(frames(i).cdata, 256);
                if i == 1
                    imwrite(A, map, filename,"gif", ...
                        "LoopCount",Inf, ...
                        "DelayTime",1/fps);
                else
                    imwrite(A, map,filename,"gif", ...
                        "WriteMode","append", ...
                        "DelayTime",1/fps);
                end
            end
        end


        function inds = get_xyIndicies(pattern, mtexHREBSD)
            inds = [0,0];
            xstep = mtexHREBSD.scan.xStep;
            ystep = mtexHREBSD.scan.yStep;
            inds(1) = mtexHREBSD.findIndexXY( ...
                [pattern.xy(1) + xstep, pattern.xy(2)] ...
            );
            inds(2) = mtexHREBSD.findIndexXY( ...
                [pattern.xy(1), pattern.xy(2) + ystep] ...
            );
        end


        function alpha = partialcurl(betaderiv1,betaderiv2)
            alpha = zeros(6,size(betaderiv1,3),size(betaderiv1,4));
            alpha(1,:,:) = betaderiv2(1,1,:,:) - betaderiv1(1,2,:,:); %13
            alpha(2,:,:) = betaderiv2(2,1,:,:) - betaderiv1(2,2,:,:); %23
            alpha(3,:,:) = betaderiv2(3,1,:,:) - betaderiv1(3,2,:,:); %33
            alpha(4,:,:) = betaderiv1(1,3,:,:); %12
            alpha(5,:,:) = -betaderiv2(2,3,:,:); %21
            alpha(6,:,:) = -betaderiv2(1,3,:,:) - betaderiv1(2,3,:,:); %11 - 22
        end


        function q = get_subpixShiftImage(ccMap)
            q_p = get_subpixShift(ccMap);
            q = [q_p(1); q_p(2)];
        end


        function q_p = get_shift(testPC_p, refPC_p, r, F)
            Fr = F*r;
            deltaP = testPC_p - refPC_p;
            q_p = deltaP - r + Fr*(-testPC_p(3))/Fr(3);
        end


        function qstar_p = reverse_PC_offset(testPC_p, refPC_p, q, r)
            deltaP = testPC_p - refPC_p;
            qstar_p = (q + r - deltaP)*refPC_p(3)/testPC_p(3) - r;
        end
    end
end