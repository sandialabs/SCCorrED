function updatedOutput = mtexHREBSD_hybrid(input, output, n0, minSize)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    updatedOutput = output;
    grainIds = input.grains.id;
%     grainIds = [578, 773, 733, 724, 508, 401, 641];
    D = parallel.pool.DataQueue;
    afterEach(D, @nUpdateWaitbar)
    wb = waitbar(0, 'Starting hybrid analysis...');
    nUpdateWaitbar(length(grainIds), wb);
    tic
    for i = 1:length(grainIds)
        grain = input.get_grain(grainIds(i));
        checkPix = length(grain.ebsd) > minSize;%checkNumGrainPixels(input, grainIds(i), n0);
        if checkPix
%             disp(num2str(i)+"/"+num2str(length(grainIds)))
%             grain = input.get_grain(grainIds(i))
            idsInGrain = grain.ebsd.id;
            ids = randsample(idsInGrain, n0, false);
%             ids = grain.get_multipleGrainRefs(n0);
            ids(1) = input.refIds(grainIds(i));
            calcFs = gethybridCalcFGrain(input, ids);
            Fbar = calcFBar(calcFs);
            deltaFbar = calcDeltaFBar(output, ids);
            F1 = deltaFbar\Fbar;
            updatedOutput = updateGrain(updatedOutput, idsInGrain, F1);
        end
        send(D,1)
    end
    toc
end


function calcFRefs = gethybridCalcFGrain(input, ids)
    calcFRefs = cell(size(ids)); 
    for i = 1:length(calcFRefs)
%         disp(num2str(i)+"/"+num2str(length(calcFRefs)))
        pat = input.get_pattern(ids(i));
        calcFRefs{i} = dynamicSimulatedCalcF(pat, input);
    end
end


function Fbar = calcFBar(calcFRefs)
    Rs = zeros([3,3,length(calcFRefs)]);
    Us = zeros([3,3,length(calcFRefs)]);
    for i = 1:length(calcFRefs)
        [Rs(:,:,i), Us(:,:,i)] = poldec(calcFRefs{i}.F);
    end
    rots = rotation.byMatrix(Rs);
    quats = quaternion(rots);
    Rbar = matrix(mean(quats));
    Uhat = pagemtimes(Rbar', pagemtimes(Rs, Us));
    Fbar = Rbar*mean(Uhat,3);
end


function deltaFbar = calcDeltaFBar(output, ids)
    Rs = zeros([3,3,length(ids)]);
    Us = zeros([3,3,length(ids)]);
    for i = 1:length(ids)
        [Rs(:,:,i), Us(:,:,i)] = poldec(output.F(:,:,ids(i)));
    end
    rots = rotation.byMatrix(Rs);
    quats = quaternion(rots);
    Rbar = matrix(mean(quats));
    Uhat = pagemtimes(Rbar', pagemtimes(Rs, Us));
    deltaFbar = Rbar*mean(Uhat,3);
end


function updatedOutput = updateGrain(output, ids, F1)
    FGrain = output.F(:,:,ids);
    FGrainUpdate = pagemtimes(FGrain, F1);
    updatedOutput = output;
    updatedOutput.F(:,:,ids) = FGrainUpdate;
    updatedOutput.beta(:,:,ids) = FGrainUpdate - eye(3);
end