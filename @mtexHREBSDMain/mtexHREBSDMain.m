classdef mtexHREBSDMain < dynamicprops
    % mtexHREBSDMain performs HREBSD for an EBSD scan
    % 
    %   Uses calcf objects to calculate the lattice distortion beta and its
    %   dertivatives using calcf objects. Is capable of asychronous
    %   parallel processing using the parallel processing toolbox, or a
    %   single core without.
    %
    % Syntax
    %   obj = mtexHREBSDMain(input)
    %   obj = mtexHREBSDMain(__,Name,Value)
    %
    % Input
    %   input - @mtexHREBSD
    %
    % Name-Value Arguments
    %   calcFStrain - defines distortion calculations
    %       [] (default)|@mtexHREBSD_calcf2 
    %       If empty, skips the lattice distortion calculations.
    %   calcFAlpha - defines distortion derivative calculations
    %       [] (default)|@mtexHREBSD_calcf2
    %       If empty, skips the distortion derivative calculations.
    %   Subset - subset of data to run
    %       [] (default)|vector of @EBSD ids
    %       Ids to perform analysis at, if empty analysis will occur at
    %       each pixel of scan.
    %   ParallelScheme - determines parallelization
    %       0 (default)|2
    %       if 0 then a single process is used, if 2 is given then
    %       asychronous parallelization is used.
    %   numCores - If a parallel pool is not active, number of cores to 
    %   intialize parallel pool.
    %       1 (default)|int
    %       Skipped if a parallel pool is already active, requires
    %       Parallelscheme not be 0.
    %   ChunkSize - *WIP: KEEP AS DEFAULT (1)!* number of calculations to 
    %   send to each worker
    %       1 (default)|int
    %       Each process will require 4xChunkSize patterns allocated to
    %       memory, a 1024x1024 pattern require ~8.5MB. Total memory
    %       allocation for patterns is numCores x ChunkSize x pattern
    %       allocation. For a chunk size of 1 using 2 processes and
    %       1024x1024 patterns memory = 2x1x8.5 = 17MB. 
    %   Verbose - displays information regarding parallelization
    %       0 (default)|1
    %       Displays the active future array to the command windowafter a 
    %       given chunk is transfered to a worker. If 0, does not display 
    %       to command windows
    %   Autosave - file location to save progress to *CURRENTLY FOR
    %   PARALLEL ONLY*
    %       []|filename
    %       Periodically saves the local output variable to the given file.
    %       If the file already exists analysis will start from where the
    %       saved output ended.
    %   AutosaveInterval - how many chunks analyzed between saves
    %       1000 (default)|double
    %
    % Properties
    %   options - structure containing the Name-Value Arguments
    %   errorLog - Cell array of the MExceptions encountered at a give 
    %       chunk by parfeval during analysis.
    %   elapsedTime - time in seconds to complete anaylsis
    %       beta (optional) - tensor components of lattice distortion as a 
    %       list of 3x3 arrays, only available if calcFStrain is given.
    %   g (optional) - sample-to-crystal rotation matricies as a list of 
    %       3x3 arrays, only available if calcFStrain is given.
    %   dbetadx (optional) - tensor components of the lattice distortion
    %       derivative with respect to x as a list of 3x3 arrays, only
    %       available if calcFAlpha is given.
    %   dbetadx (optional) - tensor components of the lattice distortion
    %       derivative with respect to y as a list of 3x3 arrays, only
    %       available if calcFAlpha is given.

    properties
        options
        errorLog
        elapsedTime (1,1) double
    end


    properties(Hidden)
        doStrain (1,1) double = 0
        doAlpha (1,1) double = 0
    end

    methods
        function obj = mtexHREBSDMain(input, options)
            arguments
                input mtexHREBSD
                options.calcFStrain {...
                    mustBeA(...
                        options.calcFStrain,...
                        ["mtexHREBSD_calcF2", "double"]...
                    )} = []
                options.calcFAlpha {...
                    mustBeA(...
                        options.calcFAlpha,...
                        ["mtexHREBSD_calcF2", "double"]...
                    )} = []
                options.ParallelScheme double {...
                    mustBeMember(options.ParallelScheme, [1,2])...
                    } = []
                options.numCores (1,1) double = 1
                options.ChunkSize (1,1) double = 1
                options.Verbose {mustBeMember(options.Verbose, [0,1])} = 0
                options.Subset = []
                options.Autosave {...
                    mustBeA(...
                        options.Autosave,...
                        ["string", "char", "double"]...
                    )} = []
                options.AutosaveInterval (1,1) double = 1000
                options.QueueSizeMult (1,1) double = 2
            end
            obj.initializeBeta(input, options.calcFStrain)
            obj.initializeAlpha(input, options.calcFAlpha)
            if ~isempty(options.Subset)
                obj.checkAddProp('Subset');
                obj.Subset = options.Subset;
            end
            obj.options = options;
            obj.errorLog = {};
            obj.main(input, options);
        end


        function obj = main(obj, input, options)
            if isempty(options.ParallelScheme) || options.numCores == 1
                obj.mainSingle(input, options);
            elseif options.ParallelScheme == 1
                obj.mainParGhosts2(input, options);
            elseif options.ParallelScheme == 2
                obj.mainParGhosts(input, options);
            end
            obj.elapsedTime = toc;
        end


        function readOutput(obj, out)
            ids = out(:,1);
            nonzero = ids > 0;
            ids = ids(nonzero);
            if isprop(obj, 'beta')
                obj.beta(:,:,ids) = reshape( ...
                    out(nonzero,2:10)', [3,3,length(ids)] ...
                    );
                obj.F = obj.beta + eye(3);
                obj.g(:,:,ids) = reshape( ...
                    out(nonzero,11:19)', [3,3,length(ids)] ...
                    );
                obj.PCInitial(:,ids) = out(nonzero, 38:40)';
                if isprop(obj, "PCCalibrated")
                    obj.PCCalibrated(:,ids) = out(nonzero, 47:49)';
                end
                obj.SSE(ids) = out(nonzero, 41);
                obj.R2(ids) = out(nonzero , 42);
            end
            if isprop(obj, 'dbetadx')
                obj.dbetadx(:,:,ids) = reshape( ...
                    out(nonzero,20:28)', [3,3,length(ids)] ...
                    );
                obj.dbetady(:,:,ids) = reshape( ...
                    out(nonzero,29:37)', [3,3,length(ids)] ...
                    );
                obj.dxSSE(ids) = out(nonzero, 43);
                obj.dxR2(ids) = out(nonzero, 44);
                obj.dySSE(ids) = out(nonzero, 45);
                obj.dyR2(ids) = out(nonzero, 46);
            end 
        end


        function initializeAlpha(obj, input, calcFAlpha)
            if ~isempty(calcFAlpha)
                obj.doAlpha = 1;
                checkAddProp(obj, 'dbetadx')
                checkAddProp(obj, 'dbetady')
                checkAddProp(obj, 'dxSSE')
                checkAddProp(obj, 'dxR2')
                checkAddProp(obj, 'dySSE')
                checkAddProp(obj, 'dyR2')
                obj.dbetadx = zeros(3,3,length(input.ebsd));
                obj.dbetady = zeros(3,3,length(input.ebsd));
            end
        end


        function initializeBeta(obj, input, calcFStrain)
            if ~isempty(calcFStrain)
                obj.doStrain = 1;
                checkAddProp(obj, 'beta')
                checkAddProp(obj, 'F')
                checkAddProp(obj, 'g')
                checkAddProp(obj, 'PCInitial')
                checkAddProp(obj, 'SSE')
                checkAddProp(obj, 'R2')
                obj.beta = zeros(3,3,length(input.ebsd));
                obj.F = zeros(3,3,length(input.ebsd));
                obj.g = zeros(3,3,length(input.ebsd));
                obj.PCInitial = zeros(3, length(input.ebsd));
                if calcFStrain.AnalysisType == 2
                    checkAddProp(obj, 'PCCalibrated')
                    obj.PCCalibrated = zeros(3, length(input.ebsd));
                end
            end
        end

        
        function checkAddProp(obj, propname)
            if ~isprop(obj, propname)
                addprop(obj, propname);
            end
        end
    end


    methods(Static)
        function ebsd = ebsdSubset(input, subset)
            if ~isempty(subset)
                ebsd = input.ebsd(subset);
            else
                ebsd = input.ebsd;
            end
        end
    end
end