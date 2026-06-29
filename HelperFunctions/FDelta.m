function [refPat, newF, pcs] = FDelta(testPat, mtexHREBSD, options)
%FDelta FDelta method for mtexHREBSD
%   Detailed explanation goes here
%     p = create_inputParser;
%     parse(p, varargin{:});
%     options = p.Results;
    if options.TroubleshootingPlots
        ts = troubleshootingFunctions;
    end
    mtexHREBSD.analysis.assumptions = 'trace=0';
    ft = mtexHREBSD.ft;
    C = testPat.C;
    CM = C.M;
    Qps = ft.Qp2s;
    maxIter = mtexHREBSD.analysis.iterationLimit;
    pcs = zeros(maxIter,3);
    SSE = zeros(maxIter,1);
    R2 = zeros(maxIter,1);
    refPat = testPat.get_simulatedPattern;
%     calcF_array = cell(1,maxIter);
    for i = 1:maxIter
        pcs(i,:) = refPat.patternCenter;
        calcF = classicCalcF(refPat, testPat, mtexHREBSD);
        SSE(i) = calcF.fit.metrics.SSE;
        R2(i) = calcF.fit.metrics.SSE;
        if options.TroubleshootingPlots
            qs = mtexHREBSD.ft.phosphorFrame2imageVec(calcF.qs);
%             ts.plot_shift(refPat, mtexHREBSD, 5*qs);
            qFit = mtexHREBSD.ft.phosphorFrame2imageVec(calcF.qFit);            
            ts.plot_twoShifts(refPat, mtexHREBSD, qs, qFit)
        end
        pc_iter = refPat.patternCenter;
        [newF, delta, res] = ResolveFandDelta(calcF.F, refPat.g, Qps, CM);
        refPat = update_refPat(refPat, newF, delta, mtexHREBSD);
        testPat = update_testPat(refPat, testPat);
%         calcF_array{i} = calcF;
        converged = check_convergence(pc_iter, refPat.patternCenter, options.Tolerance);
        if converged
            pcs = pcs(1:i,:);
            SSE = SSE(1:i);
            R2 = R2(1:i);
            break
        end
        if options.Verbose
            disp("iteration: " + num2str(i-1))
            disp(pc_iter)
            disp(newF)
            disp(delta)
            disp(res)
        end
    end
%     if ~converged
%         refPat.patternCenter = [0, 0, 0];
%     end
    if options.ConvergencePlots 
        plot_patternCenters(pcs)
    end
end


function check = check_convergence(pc1, pc2, tolerance)
    diff = abs(pc2 - pc1)./pc2;
    if all(diff < tolerance)
        check = 1;
    else
        check = 0;
    end
end


function plot_fitMetrics(SSE, R2)
    iterations = 1:size(SSE, 1);
%     titles = ["x*", "y*", "z*"];
    figure
    yyaxis left
    plot(iterations, SSE, '.')
    yyaxis right
    plot(iterations, R2, '.')
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


function newTestPat = update_testPat(refPat, testPat)
    newTestPat = testPat;
    newTestPat.patternCenter = refPat.patternCenter;
    newTestPat.rotations = refPat.rotations;
    newTestPat.g = refPat.g;
end


function newRefPat = update_refPat(refPat, F, delta, mtexHREBSD)
    copyPat = refPat;
    Qps = mtexHREBSD.ft.Qp2s;
    g = refPat.g;
    [rr, ~] = poldec(g*Qps*F*Qps'*g');
    copyPat.patternCenter = [...
        refPat.patternCenter(1) - delta(1)*refPat.patternCenter(3),...
        refPat.patternCenter(2) + delta(2)*refPat.patternCenter(3),...
        refPat.patternCenter(3) + delta(3)*refPat.patternCenter(3)];
    copyPat.g = rr'*refPat.g;
    [copyPat.rotations(1),copyPat.rotations(2),copyPat.rotations(3)] = gmat2euler(copyPat.g);
    newRefPat = copyPat.get_simulatedPattern;
end


function [F_new,Delta,res] = ResolveFandDelta(F, g, Qps, C_crystal)
%     betaM = (Qps')*(g')*F*g*Qps - eye(3);
    betaM = F - eye(3);
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


function p = create_inputParser
    p = inputParser;
    addParameter(p, 'Verbose', 0)
    addParameter(p, 'ConvergencePlots', 0)
    addParameter(p, 'TroubleshootingPlots', 0)
    addParameter(p, 'Tolerance', 1E-5)
end