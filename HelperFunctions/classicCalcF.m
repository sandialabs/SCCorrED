classdef classicCalcF < mtexHREBSD_calcF
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        shifts 
        q_cc_vec
%         q_p_vec
        qFit
        testPC_p
        refPC_p
        deltaP_p
        rRef_p
%         rois
%         roisShifted
%         shifts_p
    end


    methods
        function obj = classicCalcF(refPattern, testPattern, mtexHREBSD, varargin)
            %UNTITLED13 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@mtexHREBSD_calcF(mtexHREBSD);
            options = obj.parseOptions(varargin{:});
            obj.doRobustFit = options.doRobustFit;
            obj.C = testPattern.C.M*1E9;
            refPattern.image = double(refPattern.image);
            testPattern.image = double(testPattern.image);
%             testPattern = obj.update_testPatRotation(refPattern, testPattern, mtexHREBSD);
            obj.rRef_p = zeros(3, obj.numRois);
%             obj.F_guess = obj.ft.Qp2s'*testPattern.g'*refPattern.g*obj.ft.Qp2s;
            obj.F_guess = eye(3);
            obj.refPC_p = obj.ft.pc2phosphorFrame(refPattern.patternCenter);
            obj.testPC_p = obj.ft.pc2phosphorFrame(testPattern.patternCenter);
            obj.deltaP_p = obj.testPC_p - obj.refPC_p;
            obj = obj.get_rRef_p;
            obj.scanIndex = testPattern.scanIndex;
            obj.F_guess = options.F_guess;
            obj = obj.performAnalysis(refPattern, testPattern, mtexHREBSD);
            obj.g = refPattern.g;
            [obj.beta, obj.fit, obj.qFit, obj.A, obj.b] = obj.calc_betaDetector(mtexHREBSD, obj.rs, obj.qs);
            g = refPattern.g;
%             obj.F = obj.calc_deviatoricVR;
            obj.F = obj.beta + eye(3);
            [r,~] = poldec(g*obj.ft.Qp2s*obj.F*obj.ft.Qp2s'*g');
            obj.g = r'*refPattern.g;
        end

        
        function obj = performAnalysis(obj, refPat, testPat, mtexHREBSD)
%             obj.rois = cell(length(obj.numRois), 3);
%             obj.roisShifted = cell(length(obj.numRois), 3);
            for i = 1:obj.numRois
                ri_p = obj.rRef_p(:,i);
                shiftRoi_p = obj.get_shift(ri_p, obj.F_guess);
%                 obj.shifts_p(i,:) = shiftRoi_p;
%                 shiftRoi_p = (obj.refPC_p + obj.deltaP_p + ri_p*obj.testPC_p(3)/obj.refPC_p(3))- obj.roiCenters_p(:,i);
                shiftRoi = (obj.ft.phosphorFrame2imageVec(shiftRoi_p));
                obj.shifts(i,:) = shiftRoi;
                [refRoi, rrange, crange] = mtexHREBSD.roi.get_squareRoi(refPat, i);
%                 testRoi = mtexHREBSD.roi.get_squareRoi(testPat, i);
                [testRoi, srrange, scrange] = mtexHREBSD.roi.get_squareRoiShifted( ...
                    testPat, i, shiftRoi ...
                    );
                ccMap = get_crossCorrelation(refRoi, testRoi, ...
                    mtexHREBSD.roi, 'Gradient', mtexHREBSD.analysis.Gradient);
                obj.peakHeight = max(ccMap(:));
%                 q2 = get_subpixShift(ccMap);
%                 [~,featureShift] = max(ccMap(:));
%                 [row, col] = ind2sub(size(ccMap), featureShift);
                q = obj.get_subpixShiftImage(real(ccMap));
                obj.q_cc_vec(i,:) = q;
                q = q + [mean(scrange) - mean(crange); mean(srrange) - mean(rrange)];
%                 obj.q_cc_vec(i,:) = q;
                q_p = obj.ft.imageVec2phosphorFrame(q);
%                 obj.q_p_vec(:,i) = q_p/norm(ri_p);
                qstar_p = obj.reverse_PC_offset(q_p, ri_p);
                obj.qs(:,i) = qstar_p/norm(ri_p);
                obj.rs(:,i) = ri_p/norm(ri_p);
            end
        end


        function q_p = get_shift(obj, r, F)
            Fr = F*r;
            deltaP = obj.deltaP_p;
            q_p = deltaP - r + Fr*(-obj.testPC_p(3))/Fr(3);
        end


        function qstar_p = reverse_PC_offset(obj, q, r)
%             q_p = [q_p(2); q_p(1); 0];
            deltaP = obj.deltaP_p;
%             deltaP(1) = -deltaP(1);
            qstar_p = (q + r - deltaP)*obj.refPC_p(3)/obj.testPC_p(3) - r;
        end

        
        function obj = get_rRef_p(obj)
            obj.rRef_p = obj.roiCenters_p - obj.refPC_p;
        end


        function options = get_options(obj, varargin)
            p = obj.create_inputParser;
            parse(p, varargin{:});
            options = p.Results;
        end
    end


    methods(Static)
        function q = get_subpixShiftImage(ccMap)
            q_p = get_subpixShift(ccMap);
%             q = [q_p(2); q_p(1)];
            q = [q_p(1); q_p(2)];
        end


        function q_p = offset_estimate_simple(r, refPC, testPC)
            deltaP = testPC - refPC;
            q_p = deltaP + r*deltaP(3)/refPC(3);
        end
        

        function options = parseOptions(varargin)
            p = inputParser;
            addParameter(p, 'F_guess', eye(3))
            addParameter(p, 'doRobustFit', 0)
            parse(p, varargin{:});
            options = p.Results;
        end
    end
end