function gridIndicies = get_grid(ebsd, spacing, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    options = parse_inputs(varargin{:});
    spatialLogic = get_spatialLogic(ebsd, spacing, options);
    phaseLogic = get_phaseLogic(ebsd, options);
    totalLogic = logical(spatialLogic.*phaseLogic);
    gridIndicies = ebsd.id(totalLogic);
end


function spatialLogic = get_spatialLogic(ebsd, spacing, options)
    x = ebsd.prop.x;
    y = ebsd.prop.y;
    if isempty(options.xOffset)
        xSteps = spacing:spacing:max(x,[],'all')-spacing;
    else
        xSteps = options.xOffset:spacing:max(x,[],'all')-options.xOffset;
    end
    if isempty(options.yOffset)
        ySteps = spacing:spacing:max(y,[],'all')-spacing;
    else
        ySteps = options.yOffset:spacing:max(y,[],'all')-options.yOffset;
    end
    delta = 1e-6;
    xlogic = arrayfun(@(x) min(abs(x - xSteps)) < delta, ebsd.prop.x);
    ylogic = arrayfun(@(y) min(abs(y - ySteps)) < delta, ebsd.prop.y);
    spatialLogic = xlogic.*ylogic;
end


function phaseLogic = get_phaseLogic(ebsd, options)
    if isempty(options.phase)
        phaseLogic = ones(size(ebsd));
    else
        phaseMapId = check_mineralList(ebsd, options.phase);
        phaseLogic = ebsd.phase == phaseMapId;
    end
end



function options = parse_inputs(varargin)
    p = inputParser;
    addParameter(p,'phase', [])
    addParameter(p, 'xOffset', [])
    addParameter(p, 'yOffset', [])
    addParameter(p, 'grains', [])
    parse(p, varargin{:})
    options = p.Results;
end

% function spatialLogic = get_spatialLogic(ebsdGrid, spacing, options)
%     x = ebsdGrid.prop.x;
%     y = ebsdGrid.prop.y;
%     if isempty(options.xOffset)
%         xSteps = spacing:spacing:max(x,[],'all')-spacing;
%     else
%         xSteps = options.xOffset:spacing:max(x,[],'all')-options.xOffset;
%     end
%     if isempty(options.yOffset)
%         ySteps = spacing:spacing:max(y,[],'all')-spacing;
%     else
%         ySteps = options.yOffset:spacing:max(y,[],'all')-options.yOffset;
%     end
%     xLogic = zeros(size(x));
%     yLogic = zeros(size(x));
%     for i = 1:length(xSteps)
%         xLogicIter = abs(x - xSteps(i)) < 1E-10;
%         xLogic = xLogic + xLogicIter;
%     end
%     
%     for i = 1:length(ySteps)
%         yLogicIter = abs(y - ySteps(i)) < 1E-10;
%         yLogic = yLogic + yLogicIter;
%     end
%     spatialLogic = xLogic.*yLogic;
% end


% function phaseLogic = get_phaseLogic(ebsdGrid, options)
%     if isempty(options.phase)
%         phaseLogic = ones(size(ebsdGrid));
%     else
%         phaseMapId = check_mineralList(ebsdGrid, options.phase);
%         phaseLogic = ebsdGrid.phase == phaseMapId;
%     end
% end


function phaseMapId = check_mineralList(ebsd, phase)
    for i = 1:length(ebsd.mineralList)
        if strcmp(ebsd.mineralList(i), phase) 
            phaseMapId = ebsd.phaseMap(i);
            break
        end
    end
end