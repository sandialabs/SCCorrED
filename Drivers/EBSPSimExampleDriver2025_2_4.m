clc;clear;close all


workingdir = "E:\2024\Clones\G3";
savedir = fullfile(workingdir, "Analysis\Step2\DynamicNoGrad");
fname = "Step2\h5oina\Step2 Specimen 1 Site 1 Map Data 1.h5oina";
filenames = {'scanfile', fullfile(workingdir, fname);...
            'imagefile', fullfile(workingdir, fname)};
%
analysisOptions = {...
    "numCores", 11, ...
    'iterationLimit',20,...
    'numRois',48, ...
    'simulationType', 'Dynamic',...
    'Gradient', 1,...
    'assumption', 'free-surface'...
    };

% Grain options
grainOptions = {"misorientation", 5, "minSize", 50};

% Set up data object
input = mtexHREBSD('fileOptions', filenames,...
                    'analysisOptions', analysisOptions,...
                    'grainOptions', grainOptions);

% create EBSPSim object -- this is done to be able to use dot indexing when
% simulating a pattern, mostly implemented for multi-phase materials
simulator = EBSPSim(input, 'Ni'); % single phase
% simulators = EBSPSim(input, {'Ni', 'CuAl2'}); % multiple phases
% p = input.get_pattern(1);
% p1 = p.simulate(simulators.Ni); % simulate with simulator, dot index the phase
% p2 = p.simulate(simulators.(p.material)); % more generally
% p3 = p.simulate; % Will call from EMsoft/Kinematic function if no EBSPSim input
% 
% figure 
% p2.show_pattern;
% figure
% p3.show_pattern;

clc; close all
% How to use EBSPSim with mtexHREBSD_calcf2
calcf1 = input.initializeCalcF('Simulators',simulator, 'AnalysisType',2);
% mtexHREBSD_calcf2 without EBSPSim
calcf2 = input.initializeCalcF('AnalysisType',2);

grid = get_grid(input.ebsd,20);
% scatter(input.ebsd.prop.x(grid), input.ebsd.prop.y(grid), 'kx');
%
% grid(input.ebsd.grainId(grid) == 7) = [];
%
x = input.ebsd.prop.x(grid);
y = input.ebsd.prop.y(grid);
xBoundry = input.grains.boundary.x;
yBoundry = input.grains.boundary.y;
dist = zeros(size(x));
for i = 1:length(x)
    dist(i) = min(sqrt( (x(i)-xBoundry).^2 + (y(i)-yBoundry).^2 ));
end
grid(dist < 0.1) = [];
%%
clc
a = tic;
output = mtexHREBSDMain(input, ...
    'calcFStrain', calcf1,...
    'ParallelScheme',1,...
    'Subset',grid, ...
    'numCores', 11,...
    'ChunkSize',1, ...
    'QueueSizeMult',3, ...
    'Verbose', 1 ...
    );
toc(a)