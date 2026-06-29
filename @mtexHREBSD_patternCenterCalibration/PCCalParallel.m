function [patternCenter, initialPatternCenter, location] = ...
    PCCalParallel(obj, mtexHREBSD, kwargs)
%     arguments
%         obj mtexHREBSD_patternCenterCalibration
%         mtexHREBSD 
%         kwargs.Verbose {mustBeMember(kwargs.Verbose, [0,1])} = 0
%         kwargs.ForceInitialPC (1,3) double = [0,0,0];
%     end
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    disp("Starting PC Calibration...")
    calcf = mtexHREBSD.initializeCalcF( ...
        'AnalysisType', 2 ...
        );
%     angle = pi/2 - mtexHREBSD.scan.sampleTilt*pi/180 + ...
%         mtexHREBSD.scan.cameraElevation*pi/180;
    x = mtexHREBSD.ebsd.prop.x(obj.gridIndicies);
    y = mtexHREBSD.ebsd.prop.y(obj.gridIndicies);
%     if strcmp(mtexHREBSD.patternExt, ".h5oina")
%         initialPatternCenter = get_PCH5oina(mtexHREBSD, obj.gridIndicies);
%     else
%         initialPatternCenter = double([ ...
%             mtexHREBSD.scan.patternCenter(1)-x./mtexHREBSD.scan.phosphorSize, ...
%             mtexHREBSD.scan.patternCenter(1)-y./mtexHREBSD.scan.phosphorSize*cos(angle), ...
%             mtexHREBSD.scan.patternCenter(1)-y./mtexHREBSD.scan.phosphorSize*sin(angle)...
%         ]);
%     end
    location = [x,y];
    outPCCal = mtexHREBSDMain(mtexHREBSD,...
        'calcFStrain',calcf,...
        'Subset',obj.gridIndicies, ...
        'ParallelScheme',2, ...
        'ChunkSize',1, ...
        'Verbose', kwargs.Verbose, ...
        'numCores', mtexHREBSD.analysis.numCores ...
        );
    initialPatternCenter = outPCCal.PCInitial(:,obj.gridIndicies)';
    patternCenter = outPCCal.PCCalibrated(:,obj.gridIndicies)';
%         if ~exist('F', 'var')
%             F(i) = parfeval(pool, @Main, 1, calcf, p, i);
%         else
%             F(end+1) = parfeval(pool, @Main, 1, calcf, p, i);
%         end
%         clc
%         disp(F)
%         if length(F) >= 2*mtexHREBSD.analysis.numCores
%             wait(F(1:mtexHREBSD.analysis.numCores))
%         end
%         if sum(strcmp({F.State}, 'finished')) %&& length(F) > 2*options.numCores
%             for j = 1:sum(strcmp({F.State}, 'finished'))
%                 [idx, outIter] = fetchNext(F);
%                 patternCenter(outIter(1), :) = outIter(2:end);
%                 F(idx) = [];
%                 nUpdateWaitbar;
%             end
%         end
%     end
%     wait(F)
%     for j = 1:length(F)
%         [idx, outIter] = fetchNext(F);
%         patternCenter(outIter(1), :) = outIter(2:end);
%         F(idx) = [];
%         nUpdateWaitbar;
%     end
end

% function patternCenter = get_PCH5oina(mtexHREBSD, inds)
%     dataPath = "/1/EBSD/Data";
%     headerPath = "/1/EBSD/Header";
%     h = double(h5read(...
%         mtexHREBSD.scan.scanfile,headerPath+"/Pattern Height" ...
%         ));
%     w = double(h5read(...
%         mtexHREBSD.scan.scanfile,headerPath+"/Pattern Width" ...
%         ));
%     VHRatio = h/w;
%     xstar0 = h5read( ...
%         mtexHREBSD.scan.scanfile, ...
%         dataPath+"/Pattern Center X");
%     xstar = (xstar0(inds) - (1-VHRatio)/2)/VHRatio;
%     ystar0 = h5read( ...
%         mtexHREBSD.scan.scanfile, ...
%         dataPath+"/Pattern Center Y");
%     ystar = ystar0(inds)/VHRatio;
%     zstar0 = h5read( ...
%         mtexHREBSD.scan.scanfile, ...
%         dataPath+"/Detector Distance");
%     zstar = zstar0(inds)/VHRatio;
%     patternCenter = double([xstar, ystar, zstar]);
% end   
