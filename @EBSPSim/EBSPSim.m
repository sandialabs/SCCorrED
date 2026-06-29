classdef EBSPSim < dynamicprops
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function obj = EBSPSim(input,Phases)
            arguments
                input mtexHREBSD
                Phases {...
                    mustBeA(...
                        Phases,...
                        ["string", "char", "cell"]...
                    )} = []
            end       
            obj.initialConstruction(input, Phases)
        end


        function initialConstruction(obj, input, Phases)
            if any(strcmp(class(Phases), {'string', 'char'}))
                obj.addPhase(input, Phases)
            else
                for i = 1:length(Phases)
                    obj.addPhase(input, Phases{i})
                end
            end
        end


        function addPhase(obj, input, phaseName)
            masterPatternKey = mat2emsoft(phaseName);
            obj.checkAddProp(phaseName)
            obj.(phaseName) = PatternSimulator(input, masterPatternKey);
        end


        function checkAddProp(obj, propname)
            if ~isprop(obj, propname)
                addprop(obj, propname);
            end
        end


        function s = size(obj)
            s = size(properties(obj));
        end
    end
end