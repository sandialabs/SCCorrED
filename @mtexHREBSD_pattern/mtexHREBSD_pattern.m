classdef mtexHREBSD_pattern
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    %TODO: Move SimData property to mtexHREBSD_calcF2, require simulation
    % analysis type to initialize

    properties
        image
        rotations
        g
        patternCenter
        scanIndex
        refIndex
        grainId
        material
        C
        imageInd
        xy
    end


    properties(Hidden = true)
        SimData  
    end


    methods
        function obj = mtexHREBSD_pattern(mtexHREBSD, ind, imageInd)
            %mtexHREBSD_pattern Construct an instance of this class
            %   Detailed explanation goes here
            obj.imageInd = imageInd;
            obj.scanIndex = ind;
            obj.grainId = mtexHREBSD.ebsd.grainId(obj.scanIndex);
            obj.refIndex = mtexHREBSD.refIds(obj.grainId);
            obj.xy = obj.get_xy(mtexHREBSD);
            obj.image = double(mtexHREBSD.patterns.getPattern( ...
                mtexHREBSD.scan.imagefile, imageInd));
%             obj.image = (I(1:2:end,1:2:end)+I(2:2:end,1:2:end)+I(1:2:end,2:2:end)+I(2:2:end,2:2:end))/4;
%             obj.rotations = [mtexHREBSD.phi1(mtexHREBSD.ebsd.id == ind), ...
%                             mtexHREBSD.Phi(mtexHREBSD.ebsd.id == ind), ...
%                             mtexHREBSD.phi2(mtexHREBSD.ebsd.id == ind)];
%             obj.g = obj.get_g_matrix(mtexHREBSD.phi1(mtexHREBSD.ebsd.id == ind), ...
%                                      mtexHREBSD.Phi(mtexHREBSD.ebsd.id == ind), ...
%                                      mtexHREBSD.phi2(mtexHREBSD.ebsd.id == ind));
%             obj.rotations = [mtexHREBSD.ebsd.orientations(ind).phi1, ...
%                              mtexHREBSD.ebsd.orientations(ind).Phi, ...
%                              mtexHREBSD.ebsd.orientations(ind).phi2];
%             obj.g = obj.get_g_matrix(mtexHREBSD.ebsd.orientations(ind).phi1, ...
%                                      mtexHREBSD.ebsd.orientations(ind).Phi, ...
%                                      mtexHREBSD.ebsd.orientations(ind).phi2);
            obj.rotations = obj.getRotation(mtexHREBSD, ind);
            obj.g = obj.getg(obj.rotations);
%             obj.g = obj.get_g_matrix(mtexHREBSD.ebsd.rotations(ind).phi1, ...
%                                      mtexHREBSD.ebsd.rotations(ind).Phi, ...
%                                      mtexHREBSD.ebsd.rotations(ind).phi2);
            obj.patternCenter = double(obj.get_patternCenter(mtexHREBSD));
            [obj.material, obj.C] = obj.get_material(mtexHREBSD, ind);
            obj = obj.setSimData(mtexHREBSD);
        end


        function obj = setSimData(obj, mtexHREBSD)
            obj.SimData = struct();
            obj.SimData.simulationType = mtexHREBSD.analysis.simulationType;
            obj.SimData.Gradient = mtexHREBSD.analysis.Gradient;
            obj.SimData.iterLimit = mtexHREBSD.analysis.iterationLimit;
            obj.SimData.pixelSize = mtexHREBSD.scan.pixelSize;
            obj.SimData.cameraElevation = mtexHREBSD.scan.cameraElevation;
            if isprop(mtexHREBSD.scan,"detectorOrientation")
                obj.SimData.detectorOrientation = mtexHREBSD.scan.detectorOrientation;
            else
                obj.SimData.detectorOrientation = [0,obj.SimData.cameraElevation+90,0]'*pi/180;
            end
            obj.SimData.delta = mtexHREBSD.scan.delta;
            obj.SimData.EMdataPath = mtexHREBSD.paths.EMdataPath;
            obj.SimData.EMsoftPath = mtexHREBSD.paths.EMsoftPath;
            obj.SimData.Filter = mtexHREBSD.patterns.filter;
            obj.SimData.pixelSize = mtexHREBSD.scan.pixelSize;
            if strcmp(mtexHREBSD.analysis.simulationType, 'Kinematic')
                obj.SimData.Av = mtexHREBSD.scan.KV*1E3;
                obj.SimData.sampletilt = mtexHREBSD.scan.sampleTilt*pi/180;
                obj.SimData.elevang = mtexHREBSD.scan.cameraElevation*pi/180;
                obj.SimData.material = mtexHREBSD.readMaterial(obj.material);
                obj.SimData.ft = mtexHREBSD.ft;
                obj.SimData.F = eye(3);
            end
            [~,~,ext] = fileparts(mtexHREBSD.scan.scanfile);
            obj.SimData.ScanType = ext;
        end
        

        function xy = get_xy(obj, mtexHREBSD)
            x = mtexHREBSD.ebsd.prop.x(obj.scanIndex);
            y = mtexHREBSD.ebsd.prop.y(obj.scanIndex);
            xy = [x,y];
        end


        function obj = construct_patternNaive(obj, mtexHREBSD, ind, imageInd)
            obj.scanIndex = ind;
            obj.image = mtexHREBSD.patterns.getPattern( ...
                mtexHREBSD.scan.imagefile, imageInd);
            obj.rotations = [mtexHREBSD.phi1(mtexHREBSD.ebsd.id == ind), ...
                            mtexHREBSD.Phi(mtexHREBSD.ebsd.id == ind), ...
                            mtexHREBSD.phi2(mtexHREBSD.ebsd.id == ind)];
            obj.g = obj.get_g_matrix(mtexHREBSD.phi1(mtexHREBSD.ebsd.id == ind), ...
                                     mtexHREBSD.Phi(mtexHREBSD.ebsd.id == ind), ...
                                     mtexHREBSD.phi2(mtexHREBSD.ebsd.id == ind));
%             obj.patternCenter = mtexHREBSD.scan.patternCenter;
            obj.patternCenter = obj.get_patternCenter(mtexHREBSD);
            obj.material = obj.get_material(mtexHREBSD);
        end


        
        function patternCenter = get_patternCenter(obj, mtexHREBSD)
            if mtexHREBSD.analysis.ForceNaivePC
                ind = mtexHREBSD.ebsd.id == obj.scanIndex;
                patternCenter = obj.get_patternCenterNaive(mtexHREBSD, ind);
            else
                if string(mtexHREBSD.patternExt) == ".tiff" 
                    patternCenter = obj.get_PCTiff(mtexHREBSD);
                elseif strcmp(mtexHREBSD.patternExt, ".h5oina")
                    patternCenter = obj.get_PCH5oina(mtexHREBSD);
                else
                    ind = mtexHREBSD.ebsd.id == obj.scanIndex;
                    switch mtexHREBSD.analysis.initialPatternCenter
                        case 'naive'
                            patternCenter = obj.get_patternCenterNaive( ...
                                mtexHREBSD, ind ...
                                );
                        otherwise
                            warning("Initial pattern center method " + ...
                                mtexHREBSD.analysis.initialPatternCenter + ...
                                " unrecognized, using naive plane fit.")
                            patternCenter = obj.get_patternCenterNaive( ...
                                mtexHREBSD, ind ...
                                );
                    end
                end
            end
        end


        function patternCenter = get_PCH5oina(obj, mtexHREBSD)
            dataPath = "/1/EBSD/Data";
            headerPath = "/1/EBSD/Header";
            h = double(h5read(...
                mtexHREBSD.scan.imagefile,headerPath+"/Pattern Height" ...
                ));
            w = double(h5read(...
                mtexHREBSD.scan.imagefile,headerPath+"/Pattern Width" ...
                ));
            VHRatio = h/w;
            xstar0 = h5read( ...
                mtexHREBSD.scan.imagefile, ...
                dataPath+"/Pattern Center X", obj.scanIndex,1);
            xstar = (xstar0 - (1-VHRatio)/2)/VHRatio;
            ystar0 = h5read( ...
                mtexHREBSD.scan.imagefile, ...
                dataPath+"/Pattern Center Y", obj.scanIndex,1);
            ystar = ystar0/VHRatio;
            zstar0 = h5read( ...
                mtexHREBSD.scan.imagefile, ...
                dataPath+"/Detector Distance", obj.scanIndex,1);
            zstar = zstar0/VHRatio;
            patternCenter = [xstar, ystar, zstar];
        end        

        


        function patternCenter = get_PCTiff(obj, mtexHREBSD)
            patternCenter = zeros(1,3);
            info = imfinfo(mtexHREBSD.patterns.imageNames{obj.imageInd});
            VHRatio = info.Height/info.Width;
            xistart = strfind(info.UnknownTags.Value,'<pattern-center-x-pu>');
            xifinish = strfind(info.UnknownTags.Value,'</pattern-center-x-pu>');
            thisx = str2double(info.UnknownTags.Value(xistart+length('<pattern-center-x-pu>'):xifinish-1));
            patternCenter(1) = (thisx - (1-VHRatio)/2)/VHRatio;
            yistart = strfind(info.UnknownTags.Value,'<pattern-center-y-pu>');
            yifinish = strfind(info.UnknownTags.Value,'</pattern-center-y-pu>');
            patternCenter(2) = str2double(info.UnknownTags.Value(yistart+length('<pattern-center-y-pu>'):yifinish-1))/VHRatio;
            zistart = strfind(info.UnknownTags.Value,'<detector-distance-pu>');
            zifinish = strfind(info.UnknownTags.Value,'</detector-distance-pu>'); 
            patternCenter(3) = str2double(info.UnknownTags.Value(zistart+length('<detector-distance-pu>'):zifinish-1))/VHRatio;
        end



        function show_pattern(obj)
%             figure
            scaledImage = uint8(rescale(obj.image,0,255));
            imshow(scaledImage);
        end
        
        function simulatedPattern = simulate(obj, patsim)
            arguments
                obj mtexHREBSD_pattern
                patsim {...
                    mustBeA(...
                        patsim,...
                        ["PatternSimulator", "double"]...
                    )} = []
            end
            simulatedPattern = obj;
            if ~isempty(patsim)
                I = patsim.get_dynamic_pattern(obj);
            elseif lower(obj.SimData.simulationType) == "dynamic"
                if ~exist('Ftensor', 'var')
                Ftensor = '1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,';
                    applyDeformation = 'n';
                else 
                    applyDeformation = 'y';
                end
                I = obj.getDynamicSimulation( ...
                    Ftensor, applyDeformation ...
                    );
                I = obj.SimData.Filter.filterImage(I);
            else
                I = obj.get_kinematicSimulated;
%                 I = obj.SimData.Filter.filterImage(I);
            end
            simulatedPattern.image = I;
        end

        function simulatedPattern = get_simulatedPattern(obj, Ftensor)
            if ~exist('Ftensor', 'var')
                Ftensor = '1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,';
                applyDeformation = 'n';
            else 
                applyDeformation = 'y';
            end
%             Ftensor
            simulatedPattern = obj;
            if lower(obj.SimData.simulationType) == "dynamic"
%                 simulatedImage = gen_emsoftPattern(mtexHREBSD, obj, Ftensor, applyDeformation);
                simulatedImage = obj.getDynamicSimulation( ...
                    Ftensor, applyDeformation ...
                    );
%                 if isprop(mtexHREBSD, 'patterns')
%                     I = mtexHREBSD.patterns.filter.filterImage(simulatedImage);
%                 else
%                 filter = patterns.ImageFilter(size(simulatedImage,1));
                I = obj.SimData.Filter.filterImage(simulatedImage);
%                 end
            else
                I = obj.get_kinematicSimulated;
            end
            simulatedPattern.image = I;%simulatedImage;
        end
    end


    methods(Static)
        function rotation = getRotation(mtexHREBSD, ind)
            rotation = [mtexHREBSD.ebsd.rotations(ind).phi1, ...
                        mtexHREBSD.ebsd.rotations(ind).Phi, ...
                        mtexHREBSD.ebsd.rotations(ind).phi2];
        end
%         function patternCenter = get_PCH5oina(mtexHREBSD)
%             dataPath = "/1/EBSD/Data";
%             xstar = h5read( ...
%                 mtexHREBSD.scan.scanfile, ...
%                 dataPath+"/Pattern Center X", 1,1);
%             ystar = h5read( ...
%                 mtexHREBSD.scan.scanfile, ...
%                 dataPath+"/Pattern Center Y", 1,1);
%             zstar = h5read( ...
%                 mtexHREBSD.scan.scanfile, ...
%                 dataPath+"/Detector Distance", 1,1);
%             patternCenter = [xstar, ystar, zstar];
%         end


        function patternCenter = get_patternCenterNaive(mtexHREBSD, ind)
%             ind = find(mtexHREBSD.ebsd.id == ind);
%             theta = (mtexHREBSD.scan.sampleTilt  mtexHREBSD.scan.cameraElevation)*pi/180;
%             mtexHREBSD.scan
            ang = pi/2 - mtexHREBSD.scan.sampleTilt*pi/180 + ...
                mtexHREBSD.scan.cameraElevation*pi/180;
            xStar = mtexHREBSD.scan.patternCenter(1) - mtexHREBSD.ebsd.prop.x(ind)/mtexHREBSD.scan.phosphorSize;

            yStar = mtexHREBSD.scan.patternCenter(2) + mtexHREBSD.ebsd.prop.y(ind)/mtexHREBSD.scan.phosphorSize*cos(pi/2 - mtexHREBSD.scan.sampleTilt*pi/180 + mtexHREBSD.scan.cameraElevation*pi/180);

            zStar = mtexHREBSD.scan.patternCenter(3) + mtexHREBSD.ebsd.prop.y(ind)/mtexHREBSD.scan.phosphorSize*sin(pi/2 - mtexHREBSD.scan.sampleTilt*pi/180 + mtexHREBSD.scan.cameraElevation*pi/180);
            
            patternCenter = [xStar, yStar, zStar];
        end


        function [material, C] = get_material(mtexHREBSD, ind)
            phaseId = mtexHREBSD.ebsd.phaseId(ind);
            material = mtexHREBSD.ebsd.mineralList{phaseId};
            C = mtexHREBSD.C{phaseId};
        end
        
        function g = getg(eulerangles)
            phi1 = eulerangles(1);
            Phi = eulerangles(2);
            phi2 = eulerangles(3);
            g = zeros(length(phi1));
            cp1 = cos(phi1);
            sp1 = sin(phi1);
            cp2 = cos(phi2);
            sp2 = sin(phi2);
            cP = cos(Phi);
            sP = sin(Phi);
            g(1,1,:)= cp1.*cp2-sp1.*sp2.*cP;
            g(1,2,:) = sp1.*cp2+cp1.*sp2.*cP;
            g(1,3,:) = sp2.*sP;
            g(2,1,:)= -cp1.*sp2-sp1.*cp2.*cP;
            g(2,2,:)= -sp1.*sp2+cp1.*cp2.*cP;
            g(2,3,:)= cp2.*sP;
            g(3,1,:)=  sp1.*sP;
            g(3,2,:)= -cp1.*sP;
            g(3,3,:)=  cP;
        end


        function g = get_g_matrix(phi1, Phi, phi2)
            % Adapted from euler2gmat from Open XY:
            %euler2gmat - creates a g-matrix according to bunge for phi1,PHI,phi2 in
            %   radians
            %Accepts vectors of phi1,PHI,and phi2 (6/09/2016 Brian Jackson)
            g = zeros(length(phi1));
            cp1 = cos(phi1);
            sp1 = sin(phi1);
            cp2 = cos(phi2);
            sp2 = sin(phi2);
            cP = cos(Phi);
            sP = sin(Phi);
            g(1,1,:)= cp1.*cp2-sp1.*sp2.*cP;
            g(1,2,:) = sp1.*cp2+cp1.*sp2.*cP;
            g(1,3,:) = sp2.*sP;
            g(2,1,:)= -cp1.*sp2-sp1.*cp2.*cP;
            g(2,2,:)= -sp1.*sp2+cp1.*cp2.*cP;
            g(2,3,:)= cp2.*sP;
            g(3,1,:)=  sp1.*sP;
            g(3,2,:)= -cp1.*sP;
            g(3,3,:)=  cP;
        end
    end
end