clc
materialName = "Silicon";
abc = [5.43,5.43,5.43];
C11 = 166.0;
C12 = 64.0;
C44 = 79.6;
E = 207; % GPa
poissonRatio = 0.31;
path = 'C:\Users\wggilli\Documents\MATLAB\mtex-5.8.1\data';
fname = fullfile(path,'tensor','Silicon.GPa');
cs = crystalSymmetry('cubic',abc,'mineral',materialName);
C = stiffnessTensor.load(fname,cs)