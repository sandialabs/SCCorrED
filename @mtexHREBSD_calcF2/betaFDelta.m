function obj = betaFDelta(obj, testPat, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    obj.checkAddProp('patternCenter');
    obj.assumption = 'trace=0';
    CM = testPat.C.M;
    Qps = obj.ft.Qp2s;
    maxIter = obj.IterationLimit;
    pcs = zeros(maxIter,3);
    b = zeros([3,3,2]);
    norms = zeros([1,maxIter]);
    if ~isempty(options.MakeGIF)
        obj.checkAddProp('frames');
        numFrames = testPat.SimData.iterLimit+2;
        obj.frames = struct('cdata',[],'colormap',[]); % Why does it make me do this matlab?
        obj.frames(numFrames) = struct('cdata',[],'colormap',[]);
        testPat.show_pattern;
        xlabel(num2str(1))
        frame = getframe(gcf);
        im{1} = frame;
        obj.frames(1) = getframe(gcf);
    end
    if isprop(obj, 'Simulators')
        refPat = testPat.simulate(obj.Simulators.(testPat.material));
    else
        refPat = testPat.simulate;
    end
    pcs(1,:) = refPat.patternCenter;
    for i = 1:maxIter
        pcs(i,:) = refPat.patternCenter;
        obj.patternCenter = refPat.patternCenter;
        obj.betaClassic(testPat, refPat, options);
        obj.getIterationData(refPat, i+1, options);
        if ~isempty(options.MakeGIF)
            obj.plot(refPat)
            xlabel(num2str(i+1))
            frame = getframe(gcf);
            im{i+1} = frame;
            obj.frames(i+1) = getframe(gcf);
        end
        [newF, delta, ~] = ResolveFandDelta(obj.beta, refPat.g, Qps, CM);
        b(:,:,2) = newF - eye(3);
        norms(i) = norm(b(:,:,2)-b(:,:,1));
        b(:,:,1) = b(:,:,2);
        refPat = updateRefPat(obj, refPat, newF, delta, obj.ft);
        testPat = updateTestPat(refPat, testPat);
        if norms(i) < options.Tolerance
            pcs = pcs(1:i,:);
            obj.beta = b(:,:,2);
            obj.patternCenter = refPat.patternCenter;
            if isprop(obj,'ImageIterations')
                obj.ImageIterations(:,:, i+1:end) = [];
                obj.qIteration(:,:,i+1:end) = [];
            end
            break
        end
    end
    obj.g = obj.getg(refPat.g);
    if ~isempty(options.MakeGIF)
        close all
        obj.makeGIF(im, options.MakeGIF)
%         movie(obj.frames,10,1)
        obj.frames(find(cellfun(@isempty,{obj.frames.cdata}))) = [];
    end
    close gcf
    if options.ConvergencePlots 
        plotConvergence(norms, i)
        plot_patternCenters(pcs)
    end
end


function plot_patternCenters(patternCenters)
    iterations = 1:size(patternCenters, 1);
    titles = ["x*", "y*", "z*"];
    figure
    tiledlayout(3,1)
    for i = 1:3
        nexttile
        plot(iterations, patternCenters(:,i))
        grid on
        ylabel(titles(i))
        xlabel("Iterations")
    end
end


function plotConvergence(norms, numiters)
    figure
    iters = 1:numiters;
    toplot = norms(1:numiters);
    plot(iters, toplot)
    xlabel('Iterations')
    ylabel('Norm of beta difference')
end



function newTestPat = updateTestPat(refPat, testPat)
    newTestPat = testPat;
    newTestPat.patternCenter = refPat.patternCenter;
    newTestPat.rotations = refPat.rotations;
    newTestPat.g = refPat.g;
end


function newRefPat = updateRefPat(obj, refPat, F, delta, ft)
    copyPat = refPat;
    Qps = ft.Qp2s;
    [rr, ~] = poldec(Qps*F*Qps');
    copyPat.patternCenter = [...
        refPat.patternCenter(1) - delta(1)*refPat.patternCenter(3),...
        refPat.patternCenter(2) + delta(2)*refPat.patternCenter(3),...
        refPat.patternCenter(3) + delta(3)*refPat.patternCenter(3)];
    copyPat.g = refPat.g*rr';
    [copyPat.rotations(1),copyPat.rotations(2),copyPat.rotations(3)] = gmat2euler(copyPat.g);
%     newRefPat = copyPat.get_simulatedPattern;
    if isprop(obj, 'Simulators')
        newRefPat = copyPat.simulate(obj.Simulators.(copyPat.material));
    else
        newRefPat = copyPat.simulate;
    end
end


function [F_new,Delta,res] = ResolveFandDelta(betaM, g, Qps, C_crystal)
%     betaM = (Qps')*(g')*betaM*g*Qps - eye(3);
%     betaM = F - eye(3);
    options = optimoptions('fsolve', 'Algorithm', ...
        'Levenberg-Marquardt', 'display', 'off',...
        'FunctionTolerance', 1e-12, 'StepTolerance', 1e-8);
    [X, fval] = fsolve(@(Y) tempfun(Y,betaM, g, Qps, C_crystal),zeros(12,1), options);
    beta = [X(1) X(2) X(3);X(4) X(5) X(6);X(7) X(8) X(9)];
%     1e6*beta
    Delta = [X(10); X(11); X(12)];
%     1e6*Delta
%     F_new = g*Qps*beta*(Qps')*(g') + eye(3);
    F_new = beta + eye(3);
    res = fval(11:13);
end


function Z = tempfun(X,betaM, g, Qps, C_crystal)
    beta = [X(1) X(2) X(3);X(4) X(5) X(6);X(7) X(8) X(9)];
    Delta = [X(10); X(11); X(12)];
    M = measuredbeta(beta,Delta);
    traction = traccalc(beta, g, Qps, C_crystal);    
    Z = zeros(13,1);
    Z(1:9) = reshape(betaM - M,9,1);
    Z(10) = trace(beta);
    Z(11:13) = traction/C_crystal(1,1,1,1);
    
%     Z(6) = Z(6)/10;
%     Z(3) = Z(3)/10;
end


function traction = traccalc(beta, g, Qps, C_crystal)
    [R, U] = poldec(Qps*beta*(Qps') + eye(3));
    V = R*U*(R');
    epsilon_sample = V - eye(3);
    gact = g*(R');%(g*R*(g'))'*g;  %Could this possibly be right?  Has it ever?
    C_sample = rotate4thorder(C_crystal,gact');
    sigma_sample = fourthbysecond(C_sample, epsilon_sample);
    traction = sigma_sample*[0;0;1];
end


function M = measuredbeta(beta,Delta)
    b11 = beta(1,1);
    b12 = beta(1,2);
    b13 = beta(1,3);
    b21 = beta(2,1);
    b22 = beta(2,2);
    b23 = beta(2,3);
    b31 = beta(3,1);
    b32 = beta(3,2);
    b33 = beta(3,3);
    e1 = Delta(1);
    e2 = Delta(2);
    e3 = Delta(3);
    H = [ e3 + b11*e3 - b31*e1,      b12*e3 - b32*e1, b13*e3 - e1 - b33*e1;
        b21*e3 - b31*e2, e3 + b22*e3 - b32*e2, b23*e3 - e2 - b33*e2;
                      0,                    0,                    0];
    wd = b11/3 + b22/3 + b33/3 + (2*e3)/3 + (b11*e3)/3 + (b22*e3)/3 - (b31*e1)/3 - (b32*e2)/3;
    M = (1/(1+wd))*(beta + H - eye(3)*wd);
end


function A_prime = rotate4thorder(A,q)
    A_prime = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            for o=1:3
                                for p=1:3
                                    A_prime(i,j,k,l) = A_prime(i,j,k,l) + A(m,n,o,p)*q(i,m)*q(j,n)*q(k,o)*q(l,p);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


function sigma = fourthbysecond(C,epsilon)
    sigma = zeros(3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    sigma(i,j) = sigma(i,j) + C(i,j,k,l)*epsilon(k,l);
                end
            end
        end
    end
end