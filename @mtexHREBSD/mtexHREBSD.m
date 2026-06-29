classdef mtexHREBSD < matlab.mixin.Copyable & dynamicprops
    %mtexHREBSD Summary of this class goes here
    %   Detailed explanation goes here

    %TODO: Grain splitting method for parallel processing. For single
    %crystal/polycrystals with less grains than number of workers. Good
    %idea to have this be the standard for polycrystals to break up larger
    %grains as well. 

    properties
        ebsd
        grains
        patterns
        refIds
        roi
        ft % frame transformations
        C
        version = "0.1.0"
        isGrain = 0;
        patternCenterOffsetMatrix
        material
        patternExt
    end

    properties(SetObservable)
        patternCenterOffset
        analysis
        scan
    end

    properties(Hidden = true)
        phi1
        Phi
        phi2
        firstImage
        paths
    end

    
    methods
        function obj = mtexHREBSD(varargin)
            %mtexHREBSD Construct an instance of this class
            %   Detailed explanation goes here
            
            % Read in data from scan files. If "fileOptions" is not a given
            % input will be asked to select scan and image files
            if nargin ~= 0
                obj.paths = obj.get_paths();
                obj.scan = get_scanSettings(obj.get_fileOptions(varargin{:}));
                [~,~,obj.patternExt] = fileparts(obj.scan.imagefile);
                obj.analysis = mtexHREBSD_analysisSettings(obj.get_analysisOptions(varargin{:}));
                obj.ebsd = obj.get_ebsd;
                grainOptions = obj.get_grainOptions(varargin{:});
                obj = obj.get_calcGrains(grainOptions);
                obj.patterns = obj.makePatternProviderMtex(obj.ebsd, obj.scan.imagefile);
                obj.firstImage = obj.patterns.getPattern(obj.scan.imagefile, 1);
                obj.roi = mtexHREBSD_roiSettings(obj.analysis, obj.firstImage);
                obj = obj.update_Settings;
                obj.refIds = obj.getRefIds;
                obj.C = obj.get_stiffnesses;
                obj.ft = frameTransformations(obj);
                obj.phi1 = obj.ebsd.rotations.phi1(:);
                obj.Phi = obj.ebsd.rotations.Phi(:);
                obj.phi2 = obj.ebsd.rotations.phi2(:);
            end
        end

        
        iniWrite = iniWrite(obj, kwargs)


        function calcf = initializeCalcF(obj, kwargs)
            arguments
                obj mtexHREBSD
                kwargs.RobustFit = 0
                kwargs.AnalysisType = 0
                kwargs.Assumption = []
                kwargs.Tolerance (1,1) double = 1E-5
%                 kwargs.numStds = []
%                 kwargs.stdThreshold = []
                kwargs.MADThreshold = [];
                kwargs.Simulators {...
                    mustBeA(...
                        kwargs.Simulators,...
                        ["EBSPSim","double"]...
                    )} = []
            end
            additionalInputs = cell([1,2*size(fields(kwargs),1)]);
            keys = fields(kwargs);
            values = struct2cell(kwargs);
            for i = 1:size(fields(kwargs),1)
                additionalInputs{2*(i-1)+1} = keys{i};
                additionalInputs{2*(i-1)+2} = values{i};
            end
            calcf = mtexHREBSD_calcF2(obj, additionalInputs{:});
        end


        function s = saveobj(obj)
            props = properties(obj);
            for i = 1:length(props)
                if ~strcmp(props{i}, 'patterns')
                    s.(props{i}) = obj.(props{i});
                end
            end
        end


%         function set.scan(obj, val)
%             obj.scan = val;
%             if ~isempty(obj.scan)
%                 addlistener( ...
%                 obj, 'scan', 'PostSet', @obj.updateft ...
%                 );
%             end
%         end


        function updateft(obj, metaProp, evenData)
            if ~isempty(obj.ft)
                obj.ft = frameTransformations(obj);
            end
        end


        function set.analysis(obj, val)
            obj.analysis = val;
            if ~isempty(obj.analysis)
                addlistener( ...
                obj, 'analysis', 'PostSet', @obj.updateRoi ...
                );
            end
        end


        function updateRoi(obj, metaProp, evenData)
            obj.roi = mtexHREBSD_roiSettings(obj.analysis, obj.firstImage);
        end

        
        function attachListener(obj)
            addlistener( ...
                obj, 'patternCenterOffset', 'PostSet', @obj.pcoChange ...
                );
        end


        function set.patternCenterOffset(obj, pco)
            if ~isempty(pco)
                obj.patternCenterOffset = pco;
                obj.attachListener;
            end
        end


        function pcoChange(obj, metaProp, eventData)
            pco = obj.patternCenterOffset;
            obj.patternCenterOffsetMatrix = [...
                -1/3*pco(3),           0,     pco(1);
                          0, -1/3*pco(3),     pco(2);
                          0,           0, 2/3*pco(3)
            ];
        end


        function ind = findIndexXY(obj, xy)
            cond1 = obj.ebsd.prop.x == xy(1);
            cond2 = obj.ebsd.prop.y == xy(2);
            ind = find(cond1.*cond2);
        end


        function set.scan(obj, val)
            obj.scan = val;
            if ~isempty(obj.ft)
                obj.ft = frameTransformations(obj);
            end
        end


        function obj = get_patternCenterOffset(obj, PCCal)
            obj.patternCenterOffset = PCCal.get_patternCenterOffset(PCCal);
            pco = obj.patternCenterOffset;
            obj.patternCenterOffsetMatrix = [...
                -1/3*pco(3),           0,     pco(1);
                          0, -1/3*pco(3),     pco(2);
                          0,           0, 2/3*pco(3)
            ];
        end

        
        function C = get_stiffnesses(obj)
            phases = obj.ebsd.CSList;
            C = cell(size(phases));
            for i = 1:length(phases)
                currPhase = phases{i};
                if class(currPhase) == "crystalSymmetry"
                    if currPhase.mineral == "Face Centered Cubic"
                        fname = fullfile(obj.paths.tensorPath, "Nickel.GPa");
                    elseif currPhase.mineral == "Body Centered Cubic"
                        fname = fullfile(obj.paths.tensorPath, "Nickel.GPa");
                    elseif currPhase.mineral == "Ni"
                        fname = fullfile(obj.paths.tensorPath, "Nickel.GPa");
                    elseif currPhase.mineral == "Iron bcc (old)"
                        fname = fullfile(obj.paths.tensorPath, "iron-alpha.GPa");
                    else
                        fname = fullfile(obj.paths.tensorPath, currPhase.mineral + ".GPa");
                    end
                    if ~exist(fname,'file')
                        error("The file " + string(fname) + ...
                            " does not exist, check stiffness tensor files.")
                    else
                        C{i} = stiffnessTensor.load(fname, currPhase);
                    end
                end
            end
        end


        function C = get_C(obj, mat)
%             mat = obj.ebsd.mineral;
%             mat = "Nickel";
%             mat = "Silicon";
            matStr = string(mat);
            for i = 1:length(obj.ebsd.CSList)
                cs = obj.ebsd.CSList{i};
                if class(cs) == "crystalSymmetry"
                    if string(cs.mineral) == matStr
                        csUse = cs;
                    end
                end
            end
%             cs = obj.ebsd.CS;
            matFile = matStr + ".GPa";
            fname = fullfile(obj.paths.tensorPath, matFile);
            C = stiffnessTensor.load(fname, csUse);
        end



        function obj = replace_phase(obj, oldPhase, cs)
            newPhaseInd = length(obj.ebsd.CSList) + 1;
            newPhaseMap = max(obj.ebsd.phaseMap) + 1;
            oldPhaseInd = strcmp(obj.ebsd.mineralList, oldPhase);
            oldPhaseMap = obj.ebsd.phaseMap(oldPhaseInd);
            obj.ebsd.CSList{1, newPhaseInd} = cs;
            obj.ebsd.phaseMap(newPhaseInd) = newPhaseMap;
            obj.ebsd.phaseId(obj.ebsd.phase == oldPhaseMap) = newPhaseInd;
            obj.C = obj.get_stiffnesses;
        end


        function refIds = getRefIds(obj)
            numGrains = length(obj.grains);
            refIds = zeros(numGrains,1);
            ids = singlePixelGrains(obj.ebsd,obj.grains);
            phases = obj.ebsd.mineralList(...
                ~strcmp(obj.ebsd.mineralList, 'notIndexed')...
                );
            for i = 1:length(phases)            
            currPhase = phases{i};
            grainsPhase = obj.grains(currPhase);
                for j = 1:length(grainsPhase)
                    currGrain = grainsPhase.id(j);
                    ebsdGrain = obj.ebsd(obj.ebsd.grainId == currGrain);
                    grain = obj.grains(currGrain);
                    if ~ids(j)
                        refIds(currGrain) = grainRefId(ebsdGrain, grain);
                    else 
                        refIds(currGrain) = ebsdGrain.grainId;
                    end
                end
            end
        end


        function refIds = get_refIds(obj)
            numGrains = length(obj.grains);
            refIds = zeros(numGrains,1);
            temp = copy(obj);

%             iq = obj.ebsd.prop.iq;
%             ci = obj.ebsd.prop.ci;
%             fit = obj.ebsd.prop.fit_prias;
%             for i = 1:numGrains
%                 inGrain = find(temp.ebsd.grainId == i);
%                 grain = temp.ebsd(temp.ebsd.grainId == i);
%                 [maxIq, iqInd] = max(iq(inGrain));
%                 [maxCi, ciInd] = max(ci(inGrain));
%                 [minFit, fitInd] = min(fit(inGrain));
% 
%                 minFitTradeOff = maxIq/iq(inGrain(fitInd)) + maxCi/ci(inGrain(fitInd));
%                 maxCiTradeOff = maxIq/iq(inGrain(ciInd)) + minFit/fit(inGrain(ciInd));
%                 maxIqTradeOff = minFit/fit(inGrain(iqInd)) + maxCi/ci(inGrain(iqInd));
% 
%                 votes = [ciInd, fitInd, iqInd];
%                 [~, voteInd] = min([maxCiTradeOff, minFitTradeOff, maxIqTradeOff]);
%                 grainRefId = votes(voteInd);
%                 globalRefId = grain.id(grainRefId);
%                 refIds(i) = globalRefId;
%             end
            phases = obj.ebsd.mineralList(~strcmp(obj.ebsd.mineralList, 'notIndexed'));
            GROD = zeros(size(obj.ebsd));
            GRODNormalized = zeros(size(obj.ebsd));
%             distFromCentroid = zeros(size(obj.ebsd));
%             distFromCentroidNormalized = zeros(size(obj.ebsd));
            for i = 1:length(phases)
                
                currPhase = phases{i};
                ebsdPhase = obj.ebsd(currPhase);
                grainsPhase = obj.grains(currPhase);
                meanOrisPhase = grainsPhase.meanOrientation;
                for j = 1:length(grainsPhase)
                    currGrain = grainsPhase.id(j);
                    ebsdGrain = ebsdPhase(ebsdPhase.grainId == currGrain);
                    currMeanOri = meanOrisPhase(j);
                    GROD(ebsdGrain.id) = angle(ebsdGrain.orientations, currMeanOri);
                    GRODNormalized(ebsdGrain.id) = (1./GROD(ebsdGrain.id) - min(1./GROD(ebsdGrain.id)))/...
                        (max(1./GROD(ebsdGrain.id)) - min(1./GROD(ebsdGrain.id)));
%                     distFromCentroid(ebsdGrain.id) = obj.get_distFromCentroid(currGrain);
%                     distFromCentroidNormalized(ebsdGrain.id) = ...
%                         (1./distFromCentroid(ebsdGrain.id) - min(1./distFromCentroid(ebsdGrain.id)))/...
%                         (max(1./distFromCentroid(ebsdGrain.id)) - min(1./distFromCentroid(ebsdGrain.id)));
                end
            end
            if isfield(obj.ebsd.prop, 'bc')
                iqNormalized = (obj.ebsd.bc - min(obj.ebsd.bc))/(max(obj.ebsd.bc) - min(obj.ebsd.bc));
            else
                iqNormalized = (obj.ebsd.iq - min(obj.ebsd.iq))/(max(obj.ebsd.iq) - min(obj.ebsd.iq));
            end
            % Maximize IQ and minimize GROD
            metric = 0.5.*iqNormalized + 0.*GRODNormalized;% + 0.2.*distFromCentroidNormalized;
            for i = 1:numGrains
%                 temp = copy(obj);
                grainInds = find(temp.ebsd.grainId == i);
                grain = temp.ebsd(temp.ebsd.grainId == i);
                if ~isempty(grainInds)
                    grainMetric = metric(grainInds);
                    [~, grainRefId] = max(grainMetric);
                    globalRefId = grain.id(grainRefId);
                    refIds(i) = globalRefId;
                end
            end

        end


        function [refIdsGrain, x, y] = get_multipleGrainRefs(obj, numRefs)
%             disp(size(obj.ebsd))
            refIdsGrain = zeros(numRefs,1);
            x = zeros(numRefs,1);
            y = zeros(numRefs,1);
            kam = obj.ebsd.KAM;
            kamNormalized = (kam - min(kam))/(max(kam) - min(kam));
            iqNormalized = (obj.ebsd.iq - min(obj.ebsd.iq))/(max(obj.ebsd.iq) - min(obj.ebsd.iq));
            metric = 0.5.*iqNormalized - 0.5.*kamNormalized;
            count = 0;
            while count < numRefs
                [~, refId] = max(metric);
                metric(refId) = 0;
                x_i = obj.ebsd.prop.x(refId);
                y_i = obj.ebsd.prop.y(refId);
                dist = ((x_i-x).^2 + (y_i-y).^2).^0.5;
                if ~any(dist < 1)
                    count = count + 1;
                    x(count) = x_i;
                    y(count) = y_i;
                    globalRefId = obj.ebsd.id(refId);
                    refIdsGrain(count) = globalRefId;
                end
            end
        end


        function obj = get_calcGrains(obj, options)
            minSizeInd = find(strcmp(string(options), "minSize"));
            if ~isempty(minSizeInd)
                minSize = options{minSizeInd + 1};
                % remove minSize from options
                options([minSizeInd,minSizeInd+1])=[];
            else
                minSize = [];
            end
            minMisInd = find(strcmp(string(options), "misorientation"));
            if ~isempty(minSizeInd)
                misorientation = options{minMisInd + 1};
                % remove minSize from options
                options([minMisInd,minMisInd+1])=[];
            else
                misorientation = 10;
            end
            ebsdLocal = obj.ebsd;
            [grainsLocal, ebsdLocal.grainId] = calcGrains( ...
                ebsdLocal('indexed'),'angle',misorientation*degree, ...
                options{:});
            ebsdOut = ebsdLocal;
            if ~isempty(minSize)
                grainsOut = merge( ...
                    grainsLocal,'maxSize', minSize, 'inclusions');
            else 
                grainsOut = grainsLocal;
            end
            gid = obj.updateGID(ebsdOut, grainsOut);
            ebsdOut.grainId = gid;
            obj.grains = grainsOut;
            % try to project rotations to their fundamental region(s)
            % (this requires grain segmentation)
            try
                ebsdOut = project2FundamentalRegion(ebsdOut);
            catch me
                warning(strcat("Projection of rotation data to their ",...
                    "fundamental region(s) failed.\nThe Mtex function ",...
                    "'project2FundamentalRegion' produced the ",...
                    "following error:\n",me.message),"");
            end
            obj.ebsd = ebsdOut;
        end


        function grainHREBSD = get_grain(obj, grainID)
%             size(inds)
            grainHREBSD = copy(obj);
            grainHREBSD.ebsd = grainHREBSD.ebsd(grainHREBSD.ebsd.grainId == grainID);
            inds = grainHREBSD.ebsd.id(grainHREBSD.ebsd.grainId == grainID);
            grainHREBSD.phi1 = grainHREBSD.phi1(inds);
            grainHREBSD.Phi = grainHREBSD.Phi(inds);
            grainHREBSD.phi2 = grainHREBSD.phi2(inds);
            grainHREBSD.refIds = grainHREBSD.refIds(grainID);
            grainHREBSD.isGrain = 1;
        end


        function pattern = get_pattern(obj, ind)
            if length(ind) > 1
                for i = 1:length(ind)
                    pattern(i) = obj.get_pattern(ind(i));
                end
            else
                scanInd = obj.ebsd.id(ind);
    %             imageIndex = find(obj.indexConverter(:,3) == scanInd, 1);
    %             imageIndex = obj.indexConverter2(scanInd);
                % adjust the index if IDs were assigned according to
                % column-major order
                if obj.isColumnMajor(obj.ebsd)
                    imageIndex = obj.indexConverter2(scanInd);
                else
                    imageIndex = scanInd;
                end
                pattern = mtexHREBSD_pattern(obj, scanInd, imageIndex);
                if ~isempty(obj.patternCenterOffset)
                    pattern.patternCenter = pattern.patternCenter + obj.patternCenterOffset;
                end
                if ~isempty(obj.patternCenterOffsetMatrix)
                    M = obj.patternCenterOffsetMatrix*1/(pattern.patternCenter(3));
                    [R,~] = poldec(M + eye(3));
                    pattern.g = pattern.g * obj.ft.Qp2s * R' * obj.ft.Qp2s';
                    [rot(1), rot(2), rot(3)] = gmat2euler(pattern.g);
                    pattern.rotations = [rot(1), rot(2), rot(3)];
                end
            end
        end


        function imageInd = indexConverter2(obj,scanInd)
            [row, col] = ind2sub([obj.scan.Ny, obj.scan.Nx], scanInd);
            imageInd = sub2ind([obj.scan.Nx, obj.scan.Ny], col, row);
        end


        function indexConverter = get_indexConverter(obj)
            % If this is problematic, could do this using modulo divisions
            % with the size of data (Nx, Ny).
            x = obj.ebsd.prop.x;
            y = obj.ebsd.prop.y;
            x(abs(x) < 1E-10) = 0;
            y(abs(y) < 1E-10) = 0;
            indexConverter = [x, y, [1:length(x)]'];
            indexConverter = sortrows(indexConverter, 2);
        end


        function obj = update_Settings(obj)
            obj.scan = obj.scan.update_delta(obj.roi);
            obj.scan = obj.scan.update_pixelSize(obj.roi);
        end


        function ebsd = get_ebsd(obj)
            [~,~,ext] = fileparts(obj.scan.scanfile);
            if ext == '.ang'
                warning('off')%, 'mtex:angConversion')
                ebsd = EBSD.load( ...
                    obj.scan.scanfile, ...
                    'convertEuler2SpatialReferenceFrame', ...
                    'setting 2' ...
                    );
%                 ebsd = EBSD.load(obj.scan.scanfile);
                warning('on')%, 'mtex:angConversion')
            elseif ext == '.ctf'
                warning('off')%, 'mtex:ctfConversion')
                ebsd = EBSD.load(obj.scan.scanfile); 
                ebsd = rotate(ebsd,rotation.byAxisAngle(yvector,180*degree),'keepXY');
                warning('on')%, 'mtex:ctfConversion')
            elseif ext == '.h5oina'
                warning('off')%, 'mtex:ctfConversion')
                ebsd = loadEBSD_h5oinaNP(obj.scan.scanfile); 
                ebsd = rotate(ebsd,rotation.byAxisAngle(yvector,180*degree),'keepXY');
                warning('on')%, 'mtex:ctfConversion')
            else
                error('mtexHREBSD currently requires .ang, .ctf, .h5oina scan files')
            end
        end


        function grainOptions = get_grainOptions(obj, varargin)
            inds = obj.parse_inputs(varargin, "grainOptions");
            if ~isempty(inds)
                grainOptions = varargin{inds(2)};
            else 
                grainOptions = {};
            end
        end


        function fileOptions = get_fileOptions(obj, varargin)
            inds = obj.parse_inputs(varargin, "fileOptions");
            if ~isempty(inds)
                fileOptions = varargin{inds(2)};
            else 
                fileOptions = {};
            end
        end


        function analysisOptions = get_analysisOptions(obj, varargin)
            inds = obj.parse_inputs(varargin, "analysisOptions");
            if ~isempty(inds)
                analysisOptions = varargin{inds(2)};
            else 
                analysisOptions = {};
            end
            if strcmp(obj.patternExt, ".tiff")
                analysisOptions{end+1} = 'initialPatternCenter';
                analysisOptions{end+1} = '.tiff';
            end

        end


        function [mu, sigma] = get_propStats(obj, ebsdProp)
            data = obj.ebsd.prop.(ebsdProp);
            mu = mean(data);
            sigma = std(data);
        end


        function plot_propHist(obj, ebsdProp)
            data = obj.ebsd.prop.(ebsdProp);
            figure;
            histogram(data, 25);
        end

        
        function [M, filename] = readMaterial(obj, Material)
            % Reads material data from text file in /Materials subfolder into a
            % structure
            obj.scan.material = Material;
            if strcmp(Material, 'notIndexed')
                M = {};
                filename = [];
                return
            end
            if strcmp(Material,'newphase')
                Material='iron-alpha';
            end
            if strcmp(Material,'Ni')
                Material='nickel';
            end
            if strcmpi(Material,'Austenite')
                Material='iron-gamma';
            end
            if strcmpi(Material,'Ferrite')
                Material='iron-alpha';
            end
            if strcmpi(Material,'aluminium')
                Material = 'aluminum';
            end
            
            filename = fullfile(obj.paths.materialPath,lower(Material)+".txt");
            
            if exist(filename,'file')
                fid = fopen(filename);
                % Import each line into a structure with field names determined by file
                while ~feof(fid)
                    tline = fgetl(fid);
                    [parname, value] = strtok(tline);
                    if ~strcmp(parname,'Material') && ~strcmp(parname,'lattice')
                        value = sscanf(value,'%f');
                    else
                        value = strtrim(value);
                    end
                    M.(parname)= value;
                end
                % Reshape hkl with correct number of columns
                if isfield(M,'hkl')
                    if strcmp(M.lattice, 'hexagonal')
                        LatticeNumber = 4;
                    else
                        LatticeNumber = 3;
                    end
                    M.hkl = reshape(M.hkl,LatticeNumber,[])';
                    if strcmp(M.lattice, 'hexagonal')
                        hkil = M.hkl;
                        hkl = zeros(size(hkil,1),3);
                        for i = 1:size(hkil,1)
                            hkl(i,1) = 3/2*hkil(i,1);
                            hkl(i,2) = sqrt(3)/2*(hkil(i,1)+2*hkil(i,2));
                            hkl(i,3) = 3/2*1/(M.c1/M.a1)*hkil(i,4);
                        end
                        M.hkl_hex = M.hkl;
                        M.hkl = hkl;
                    end
                end
                fclose(fid);
            else
                warning("Material file not found at "+filename)
                M = {};
            end
            
            %SplitDD info
            switch lower(Material)
                case 'nickel'
                    M.SplitDD = {'Ni','Ni(18ss)'};
                case 'magnesium'
                    M.SplitDD = {'Mg','Mg (a systems only)'};
                case 'copper'
                    M.SplitDD = {'Cu'};
                case 'tantalum'
                    M.SplitDD = {'Ta','Ta (with 112 planes)'};
                case 'aluminum'
                    M.SplitDD = {'Al-18ss'};
                case 'iron-alpha'
                    M.SplitDD = {'Fe'};
                case 'zirconium (alpha)'
                    M.SplitDD = {'Zr'};
            end
        end


    function distFromCentroid = get_distFromCentroid(obj, grainId)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
            grain_i = obj.grains(grainId);
            ebsd_i = obj.ebsd(obj.ebsd.grainId == grain_i.id);
            centroid_i = grain_i.centroid;
%             distFromCentroid = zeros(length(ebsd_i), 1);
            deltaX = ebsd_i.prop.x - centroid_i(1);
            deltaY = ebsd_i.prop.y - centroid_i(2);
            distFromCentroid = mean((deltaX.^2 - deltaY.^2).^0.5,1);
%             for j = 1:length(ebsd_i)
%                deltaX_j = ebsd_i.prop.x(j) - centroid_i(1);
%                deltaY_j = ebsd_i.prop.y(j) - centroid_i(2);
%                distFromCentroid(j) = mean((deltaX_j.^2 - deltaY_j.^2).^0.5);
%             end
    end

    end



    methods(Static)
        function gid = updateGID(ebsd, grains)
            gid = ebsd.grainId;
            x = ebsd.prop.x;
            y = ebsd.prop.y;
            for i = 1:length(grains.id)
                whichGrain = grains.id(i);
                inGrain = checkInside( ...
                    grains(whichGrain), [x,y], 'includeBoundary' ...
                    );
                gid(inGrain) = whichGrain;
            end
        end


        function obj = loadobj(s)
            obj = mtexHREBSD;
            obj.paths = obj.get_paths();
            props = fields(s);
            for i = 1:length(props)
                obj.(props{i}) = s.(props{i});
            end
            if ~isfile(obj.scan.imagefile)
                disp("Could not find "+obj.scan.imagefile+", input new image location");
                name=input("New file: " ,'s');
                obj.scan.imagefile = name;
            end
            obj.phi1 = obj.ebsd.rotations.phi1(:);
            obj.Phi = obj.ebsd.rotations.Phi(:);
            obj.phi2 = obj.ebsd.rotations.phi2(:);
            obj.patterns = obj.makePatternProviderMtex( ...
                    obj.ebsd, obj.scan.imagefile ...
                    );
            p = obj.get_pattern(1);
            obj.firstImage = p.image;
            obj.updateRoi;
            obj.ft = frameTransformations(obj);
            addlistener( ...
                obj, 'analysis', 'PostSet', @obj.updateRoi ...
                );
%             addlistener( ...
%                 obj, 'scan', 'PostSet', @obj.updateft ...
%                 );
        end


        function patternProvider = makePatternProviderMtex(ebsd, first_image)
            [~, ~, ext] = fileparts(first_image);
            switch ext
                case '.h5'
                    patternProvider =...
                        patterns.H5PatternProvider(first_image);
                case '.h5oina'
                    patternProvider =...
                        patterns.H5oinaPatternProvider(first_image);
                case {'.up1', '.up2'}
                    patternProvider =...
                        patterns.UPPatternProvider(first_image);
%                 case {'.ebsp'}
%                     Nx = length(unique(ebsd.prop.x));
%                     Ny = length(unique(ebsd.prop.x));
%                     patternProvider =...
%                         patterns.EBSPPatternProvider(first_image, Nx, Ny);
                case { '.jpg', '.jpeg', '.tif', '.tiff', '.bmp', '.png'}
                    x = single(ebsd.prop.x);
                    y = single(ebsd.prop.y);
                    x(abs(x) < 1E-5) = 0;
                    y(abs(y) < 1E-5) = 0;
                    X = unique(x);
                    Y = unique(y);
                    xStep = X(2) - X(1);
                    yStep = Y(2) - Y(1);
                    Nx = length(X);
                    Ny = length(Y);
                    patternProvider =...
                        patterns.ImagepatternProvider(...
                        first_image,...
                        'Square',...
                        length(ebsd),...
                        [Nx, Ny],...
                        [ebsd.prop.x(1), ebsd.prop.x(2)],...
                        [xStep, yStep]);
                    
                otherwise
                    error('PatternProvider:UnrecognizedExtention', ...
                        'Unrecognized file extention %s', ext)
            end
        end


        function out = parse_inputs(input, which_option)
            ind = find(strcmp(input, string(which_option)));
            out = [ind, ind+1];
        end


        function paths = get_paths()
            configFile = fullfile( ...
                getenv("SYSTEMDRIVE"), ...
                "Users", ...
                getenv('username'), ...
                ".config\mtexHREBSD\mtexHREBSDConfig.json" ...
                );
            fid = fopen(configFile);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            paths = jsondecode(str);
        end

        function colMajor = isColumnMajor(ebsd)
            % determine whether the Mtex EBSD object belonging to 'obj' 
            % traverses along rows or columns
            colMajor = true;
            try
                if ebsd(2).prop.x - ebsd(1).prop.x > 0
                    colMajor = false;
                end
            catch me
                warning(me.message);
                warning("Detection of EBSD data storage order " + ...
                    "failed. Column-major order is assumed.");
            end
        end

    end
end