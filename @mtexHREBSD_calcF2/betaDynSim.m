function obj = betaDynSim(obj, testPat, options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if isprop(obj, 'Simulators')
        refPat = testPat.simulate(obj.Simulators.(testPat.material));
    else
        refPat = testPat.simulate;
    end
    b = zeros([3,3,2]);
    obj.betaClassic(testPat, refPat, options);
    if ~isempty(options.MakeGIF)
        obj.checkAddProp('frames');
        numFrames = testPat.SimData.iterLimit+2;
        obj.frames = struct('cdata',[],'colormap',[]); % Why does it make me do this matlab?
        obj.frames(numFrames) = struct('cdata',[],'colormap',[]);
        testPat.show_pattern;
        xlabel(num2str(1))
        frame = getframe(gcf);
%         close gcf
        im{1} = frame;
        obj.frames(1) = frame;
        plot(obj, refPat)
        xlabel(num2str(2))
        frame = getframe(gcf);
        im{2} = frame;
        obj.frames(2) = frame;
%         close gcf
    end
    b(:,:,1) = obj.beta;
    norms = zeros([1,testPat.SimData.iterLimit]);
    converged = 0;
    for i = 1:obj.IterationLimit
        refPat = updateRefPat(obj,refPat, b(:,:,1));
        testPat = updateTestPat(refPat, testPat);
        obj.betaClassic(testPat, refPat, options);
        if ~isempty(options.MakeGIF)
            plot(obj, refPat)
            xlabel(num2str(i+2))
            frame = getframe(gcf);
            im{i+2} = frame;
            obj.frames(i+2) = frame;
        end
        b(:,:,2) = obj.beta;
        norms(i) = norm(b(:,:,2) - b(:,:,1));
        if norms(i) < options.Tolerance
            converged = 1;
            break
        else
            b(:,:,1) = b(:,:,2);
        end
    end
    obj.g = obj.getg(refPat.g);
    % remove empty frames
    if ~isempty(options.MakeGIF)
        if i + 3 <= length(obj.frames)
            obj.frames(i+3:end) = [];
        end
    end

    if ~converged
        disp( ...
            "Dynamic simulations did not converge for id " + ...
            num2str(testPat.scanIndex) ...
            )
    end
    if ~isempty(options.MakeGIF)
        obj.makeGIF(obj.frames, options.MakeGIF, options.GIFfps);
        obj.frames(find(cellfun(@isempty,{obj.frames.cdata}))) = [];
    end
    close gcf
    if options.ConvergencePlots
        plotConvergence(norms, i)
    end
end


function plotConvergence(norms, numiters)
    f = figure;
    iters = 1:numiters;
    toplot = norms(1:numiters);
    plot(iters, toplot)
    xlabel('Iterations')
    ylabel('Norm of beta difference')
end


function newTestPat = updateTestPat(refPat, testPat)
    newTestPat = testPat;
    newTestPat.rotations = refPat.rotations;
    newTestPat.g = refPat.g;
end

function newRefPat = updateRefPat(obj,refPat, beta)
    copyPat = refPat;
    % transform components in to the sample frame
    Qps = obj.ft.Qp2s;
    % compute rotation using polar decomposition
    F = eye(3) + Qps*beta*Qps';
    [rr, ~] = poldec(F);
    % update the orientation matrix (g)
    copyPat.g = refPat.g*rr';
    [copyPat.rotations(1),copyPat.rotations(2),copyPat.rotations(3)] =...
        gmat2euler(copyPat.g);
%     newRefPat = copyPat.get_simulatedPattern;
    if isprop(obj, 'Simulators')
        newRefPat = copyPat.simulate(obj.Simulators.(copyPat.material));
    else
        newRefPat = copyPat.simulate;
    end
end
