classdef (Abstract) mtexHREBSD_calcF
    %mtexHREBSD_calcF is the super class for calcF methods
    %   Detailed explanation goes here

    properties
        logMessage
        F_guess
        peakHeight
        ft 
        qs
        rs
        roiCenters
        roiCenters_p
        numRois
        beta
        F
        C
        g
        fit
        A
        b
        doRobustFit {mustBeInteger}
    end

    methods
        function obj = mtexHREBSD_calcF(mtexHREBSD)
            %mtexHREBSD_calcF Construct an instance of this class
            %   Detailed explanation goes here
            obj.ft = mtexHREBSD.ft;
            obj.fit = mtexHREBSD_fit("SSE", "R2");
%             obj.C = mtexHREBSD.C.M*1E9;
            obj.numRois = mtexHREBSD.roi.numRois;
            obj = obj.intializeArray;
            obj.roiCenters = mtexHREBSD.roi.centers;
            obj.roiCenters_p = obj.ft.imageVec2phosphorFrame(obj.roiCenters);
        end
        

        function out = calc_deviatoricVR(obj)
            [V,R] = polarDecompositionVR(obj.beta + eye(3,3));
            I = eye(3,3);
            out = (I + V - trace(V)*I./3.0)*R;
        end


        function [A,b] = construct_Ab(obj,mtexHREBSD)
            r1 = obj.rs(1,:)';
            r2 = obj.rs(2,:)';
            r3 = obj.rs(3,:)';
            q1 = obj.qs(1,:)';
            q2 = obj.qs(2,:)';
            zerovec = zeros(size(r1));
            A1 = [r1.*r3, r2.*r3, r3.*r3, zerovec, zerovec, zerovec, -(r1.*r1 + q1.*r1), -(r1.*r2 + q1.*r2), -(r1.*r3 + q1.*r3)];
            A2 = [zerovec, zerovec, zerovec, r1.*r3, r2.*r3, r3.*r3, -(r2.*r1 + q2.*r1), -(r2.*r2 + q2.*r2), -(r2.*r3 + q2.*r3)];
            b1 = q1.*r3;
            b2 = q2.*r3;
            C_rotate = obj.ft.Qp2s'*obj.g';
            C_p = rotate4thorder(obj.C, C_rotate);
            n = obj.ft.Qp2s'*[0;0;1];
            switch mtexHREBSD.analysis.assumptions
                case 'free-surface'
                    A5=[ C_p(1,1,1,1)*n(1)+C_p(1,2,1,1)*n(2)+C_p(1,3,1,1)*n(3) C_p(1,1,1,2)*n(1)+C_p(1,2,1,2)*n(2)+C_p(1,3,1,2)*n(3)    C_p(1,1,1,3)*n(1)+C_p(1,2,1,3)*n(2)+C_p(1,3,1,3)*n(3) C_p(1,1,1,2)*n(1)+C_p(1,2,1,2)*n(2)+C_p(1,3,1,2)*n(3) C_p(1,1,2,2)*n(1)+C_p(1,2,2,2)*n(2)+C_p(1,3,2,2)*n(3)  C_p(1,1,2,3)*n(1)+C_p(1,2,2,3)*n(2)+C_p(1,3,2,3)*n(3)  C_p(1,1,1,3)*n(1)+C_p(1,2,1,3)*n(2)+C_p(1,3,1,3)*n(3)      C_p(1,1,2,3)*n(1)+C_p(1,2,2,3)*n(2)+C_p(1,3,2,3)*n(3) C_p(1,1,3,3)*n(1)+C_p(1,2,3,3)*n(2)+C_p(1,3,3,3)*n(3)]/1e11;
                    A6=[ C_p(2,1,1,1)*n(1)+C_p(2,2,1,1)*n(2)+C_p(2,3,1,1)*n(3) C_p(2,1,1,2)*n(1)+C_p(2,2,1,2)*n(2)+C_p(2,3,1,2)*n(3)    C_p(2,1,1,3)*n(1)+C_p(2,2,1,3)*n(2)+C_p(2,3,1,3)*n(3) C_p(2,1,1,2)*n(1)+C_p(2,2,1,2)*n(2)+C_p(2,3,1,2)*n(3) C_p(2,1,2,2)*n(1)+C_p(2,2,2,2)*n(2)+C_p(2,3,2,2)*n(3)  C_p(2,1,2,3)*n(1)+C_p(2,2,2,3)*n(2)+C_p(2,3,2,3)*n(3)  C_p(2,1,1,3)*n(1)+C_p(2,2,1,3)*n(2)+C_p(2,3,1,3)*n(3)      C_p(2,1,2,3)*n(1)+C_p(2,2,2,3)*n(2)+C_p(2,3,2,3)*n(3) C_p(2,1,3,3)*n(1)+C_p(2,2,3,3)*n(2)+C_p(2,3,3,3)*n(3)]/1e11;
                    A7=[ C_p(3,1,1,1)*n(1)+C_p(3,2,1,1)*n(2)+C_p(3,3,1,1)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3)    C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3) C_p(3,1,2,2)*n(1)+C_p(3,2,2,2)*n(2)+C_p(3,3,2,2)*n(3)  C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3)  C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3)      C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3) C_p(3,1,3,3)*n(1)+C_p(3,2,3,3)*n(2)+C_p(3,3,3,3)*n(3)]/1e11;
                    A = [A1; A2; A5; A6; A7];
                    b = [b1; b2;  0;  0;  0];
                case 'traction-free'
                    A7=[ C_p(3,1,1,1)*n(1)+C_p(3,2,1,1)*n(2)+C_p(3,3,1,1)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3)    C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3) C_p(3,1,2,2)*n(1)+C_p(3,2,2,2)*n(2)+C_p(3,3,2,2)*n(3)  C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3)  C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3)      C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3) C_p(3,1,3,3)*n(1)+C_p(3,2,3,3)*n(2)+C_p(3,3,3,3)*n(3)]/1e11;
                    A = [A1; A2; A7];
                    b = [b1; b2; 0];
                case 'trace=0'
                    A7=[1 0 0 0 1 0 0 0 1];
                    b7 = 0;
                    A = [A1; A2; A7];
                    b = [b1; b2; b7];
            end
        end


        function [beta, fit, qFit, A4, b4] = calc_betaDetector(obj, mtexHREBSD)
            r1 = obj.rs(1,:)';
            r2 = obj.rs(2,:)';
            r3 = obj.rs(3,:)';
            q1 = obj.qs(1,:)';
            q2 = obj.qs(2,:)';
            [A4, b4] = obj.construct_Ab(mtexHREBSD);
%             Anorm = normalize(A4);
%             bnorm = normalize(b4);
            if obj.doRobustFit
                [X,~] = robustfit(A4,b4,'bisquare',4.685,'off');
            else
                X = A4\b4;
            end
%             obj.fit.stats = stats;
%             X = A4\b4
            beta = reshape(X, [3 3])';
            bFit = A4*X;
            q1Fit = bFit(1:length(q1))./r3;
            fitInd = obj.numRois*2;
            q2Fit = bFit(length(q1)+1:fitInd)./r3;
            qFit = [q1Fit, q2Fit]';
            qFit(3,:) = 0;
            fit = obj.fit.calc_metrics(b4, bFit);
            
        end

        function obj = intializeArray(obj)
            obj.rs = zeros(3, obj.numRois);
            obj.qs = zeros(3, obj.numRois);
%             obj.rRef_p = zeros(3, obj.numRois);
%             obj.rTest = zeros(3, obj.numRois);
        end


        function plot(obj, pattern, varargin)
            p = inputParser;
            addParameter(p, 'multiplier',1)
            addParameter(p, 'showFit',0)
            addParameter(p, 'numRois', obj.numRois)
            addParameter(p, 'showRois', 0)
            addParameter(p, 'roiSize', 0)
            parse(p, varargin{:});
            options = p.Results;
            m = options.multiplier;
            q_iv = obj.ft.phosphorFrame2imageVec(obj.qs);
            numSkip = round(obj.numRois/options.numRois);
            centers = obj.roiCenters(:, 2:numSkip:obj.numRois);
            qPlot = q_iv(:, 2:numSkip:obj.numRois);
            pattern.show_pattern;
            hold on
            q = quiver( ...
                centers(1,:), centers(2,:), ...
                m*qPlot(1,:), ...
                m*qPlot(2,:), ...
                'b', 'LineWidth',1, ...
                'DisplayName', 'Measured' ...
                );
            q.AutoScale = 'off';
%             if options.showFit
            qfit_iv = obj.ft.phosphorFrame2imageVec(obj.qFit);
            qfitPlot = qfit_iv(:, 2:numSkip:obj.numRois);
            qfit = quiver( ...
                centers(1,:), centers(2,:), ...
                m*qfitPlot(1,:), ...
                m*qfitPlot(2,:), ...
                'r','LineWidth',1, ...
                'DisplayName', 'Fit' ...
                );
            qfit.AutoScale = 'off';
            legend
%             end
            if options.roiSize
%                 x = [centers(1,:)-options.roiSize/2; centers(1,:)+options.roiSize/2]
%                 y = [centers(2,:)-options.roiSize/2; centers(2,:)+options.roiSize/2];
                for i = 1:length(centers)
                    x = [centers(1,i)-options.roiSize/2; centers(1,i)+options.roiSize/2];
                    y = [centers(2,i)-options.roiSize/2; centers(2,i)+options.roiSize/2];
                    plot([x(1), x(2), x(2), x(1), x(1)], [y(1), y(1), y(2), y(2), y(1)], 'Color', 'white')
                end
            end
            hold off
            axis equal
        end


%         function obj = get_rRef_p(obj)
%             obj.rRef_p = obj.roiCenters_p - repmat(obj.refPC_p,1,obj.numRois);
%         end
    end
end