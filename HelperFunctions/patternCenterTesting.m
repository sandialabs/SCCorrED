savedir = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\userScripts\Gilliland\Data\DevelopmentTesting';
load(fullfile(savedir,"test_angStep3.mat"))
load(fullfile(savedir,"test_angStep3PCCal.mat"))
% load(test_angStep3PCCal)
%%
patInd = 100000;
patNoCal = test_angStep3.get_pattern(patInd)
patCal = test_angStep3PCCal.get_pattern(patInd)
patCalFDelta = FDelta(patCal, test_angStep3PCCal);%, 'Verbose', 1, 'ConvergencePlots',1)
% calcF = dynamicSimulatedCalcF(patCalFDelta,  test_angStep3PCCal, varargin)


patNoCal.show_pattern
patCal.show_pattern
patCalFDelta.show_pattern
%%
patInd = 157258;
patCal = test_angStep3PCCal.get_pattern(patInd);
raw = test_angStep3.get_pattern(patInd);
[FDelta1, FDelta1F] = FDelta(raw, test_angStep3, 'ConvergencePlots', 1);
% FDelta2 = FDelta(FDelta1, test_angStep3);
% FDelta1CalcF = classicCalcF(FDelta1, raw, test_angStep3);
FDelta1.image = raw.image;
% [FDelta2, FDelta2F] = FDelta(FDelta1, test_angStep3);
dynSim1 = dynamicSimulatedCalcF(FDelta1, test_angStep3);
dynSim2 = dynamicSimulatedCalcF(raw, test_angStep3);
raw2 = raw;
raw2.patternCenter = FDelta1.patternCenter;
dynSim3 = dynamicSimulatedCalcF(raw2, test_angStep3);
% FDelta2CalcF = classicCalcF(FDelta2, FDelta1, test_angStep3);
dynSim4 = dynamicSimulatedCalcF(patCal, test_angStep3);
dynSim5 = dynamicSimulatedCalcF(raw, test_angStep3);

clc
% disp("FDelta 1 g = ")
% disp(FDelta1.g)
disp("FDelta 1 F = ")
disp(FDelta1F)
% 
% disp("")
% disp("FDelta 2 g = ")
% disp(FDelta2.g)
% disp("FDelta 2 F = ")
% disp(FDelta2F)

disp("")
% disp("Dyn sim 1 g = ")
% disp(dynSim1.g)
disp("Dyn sim 1 F = ")
disp(dynSim1.F)

% disp("")
% disp("Dyn sim 2 g = ")
% disp(dynSim2.g)
% disp("Dyn sim F 2 = ")
% disp(dynSim2.F)
% 
disp("")
% disp("Dyn sim 3 g = ")
% disp(dynSim3.g)
disp("Dyn sim F 3 = ")
disp(dynSim3.F)


disp("")
% disp("Dyn sim 4 g = ")
% disp(dynSim4.g)
disp("Dyn sim 4 F = ")
disp(dynSim4.F)

disp("")
% disp("Dyn sim 5 g = ")
% disp(dynSim5.g)
disp("Dyn sim 5 F = ")
disp(dynSim5.F)

% raw.show_pattern
% FDelta1.show_pattern
%%
dynSim4 = dynamicSimulatedCalcF(FDelta1, test_angStep3);
dynSim5 = dynamicSimulatedCalcF(raw2, test_angStep3);
disp("")
disp("Dyn sim 4 g = ")
disp(dynSim4.g)
disp("Dyn sim 4 F = ")
disp(dynSim4.F)

disp("")
disp("Dyn sim 5 g = ")
disp(dynSim5.g)
disp("Dyn sim 4 5 = ")
disp(dynSim5.F)